/**
 * 
 */
package pt.uminho.ceb.biosystems.merlin.gpr.rules.core;

import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Observable;
import java.util.Observer;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.kegg.KeggAPI;
import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.ncbi.NcbiAPI;
import pt.uminho.ceb.biosystems.merlin.database.connector.databaseAPI.HomologyAPI;
import pt.uminho.ceb.biosystems.merlin.database.connector.databaseAPI.ModelAPI;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.Connection;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.DatabaseAccess;
import pt.uminho.ceb.biosystems.merlin.local.alignments.core.RunSimilaritySearch;
import pt.uminho.ceb.biosystems.merlin.utilities.DatabaseProgressStatus;
import pt.uminho.ceb.biosystems.merlin.utilities.Enumerators.AlignmentScoreType;
import pt.uminho.ceb.biosystems.merlin.utilities.Enumerators.Method;
import pt.uminho.ceb.biosystems.merlin.utilities.TimeLeftProgress;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.capsules.AlignmentCapsule;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.gpr.ReactionProteinGeneAssociation;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.gpr.ReactionsGPR_CI;
import pt.uminho.ceb.biosystems.merlin.utilities.datastructures.map.MapUtils;
import pt.uminho.ceb.biosystems.merlin.utilities.io.FileUtils;

/**
 * @author ODias
 *
 */
public class IdentifyGenomeSubunits extends Observable implements Observer {

	private static final Logger logger = LoggerFactory.getLogger(IdentifyGenomeSubunits.class);

	private Map<String, List<String>> ecNumbers;
	private Map<String, AbstractSequence<?>> genome;
	private long reference_organism_id;
	private ConcurrentHashMap<String, AbstractSequence<?>> sequences;
	private ConcurrentHashMap<String, Set<String>> closestOrtholog;
	private DatabaseAccess dba;
	private double similarity_threshold;
	private Method method;
	private AtomicBoolean cancel;
	private double referenceTaxonomyThreshold;
	private boolean compareToFullGenome;
	private TimeLeftProgress progress;
	private ConcurrentLinkedQueue<AlignmentCapsule> findGapsResult;
	private long startTime;
	private String wsTaxonomyTempFolderPath;



	/**
	 * Constructor for class that identifies genome subunits.
	 * 
	 * @param ec_numbers
	 * @param genome
	 * @param reference_organism_id
	 * @param dba
	 * @param similarity_threshold
	 * @param referenceTaxonomyThreshold
	 * @param method
	 * @param compareToFullGenome
	 * @param cancel
	 */
	public IdentifyGenomeSubunits(Map<String, List<String>> ec_numbers, Map<String, AbstractSequence<?>> genome, long reference_organism_id, 
			DatabaseAccess dba, double similarity_threshold, double referenceTaxonomyThreshold, Method method, 
			boolean compareToFullGenome) {

		this.ecNumbers = ec_numbers;
		this.genome = genome;
		this.reference_organism_id = reference_organism_id;
		this.dba = dba;
		this.similarity_threshold = similarity_threshold;
		this.method = method;
		this.referenceTaxonomyThreshold = referenceTaxonomyThreshold;
		this.compareToFullGenome = compareToFullGenome;
		this.findGapsResult = new ConcurrentLinkedQueue<AlignmentCapsule>();
		startTime = GregorianCalendar.getInstance().getTimeInMillis();

	}


	/**
	 * @throws Exception
	 */
	public boolean runIdentification(boolean gapsIdentification) throws Exception {

		boolean ret = true;

		try {

			Connection conn = new Connection(this.dba);
			Statement statement = conn.createStatement();

			if(!this.ecNumbers.isEmpty()) {

				this.sequences = new ConcurrentHashMap<>();
				this.closestOrtholog = new ConcurrentHashMap<>();

				List<String> referenceTaxonomy = NcbiAPI.getReferenceTaxonomy(reference_organism_id);
				logger.info("Reference taxonomy set to {}", referenceTaxonomy);

				ConcurrentHashMap<String, Integer> ncbi_taxonomy_ids = new ConcurrentHashMap<>();
				ConcurrentHashMap<String, Integer> kegg_taxonomy_scores = new ConcurrentHashMap<>();
				kegg_taxonomy_scores.put("noOrg", 0);
				ConcurrentHashMap<String, Map<String, List<String>>> orthologsSequences = new ConcurrentHashMap<>();;
				Map<String, String> kegg_taxonomy_ids = IdentifyGenomeSubunits.getKeggTaxonomyIDs();

				Set<String> bypass =  ModelAPI.getECNumbersWithModules(conn);	
				List<String> iterator = new ArrayList<>(this.ecNumbers.keySet());
				iterator.removeAll(bypass);

				logger.info("Iterator size: {}, entries {}", iterator.size(), iterator);

				Map<String, Integer> geneIds = ModelAPI.getGeneIds(statement);
				Map<String, List<String>> sequenceIdsSet = ModelAPI.getSequenceIds(statement);

				GetClosestOrhologSequence seq = new GetClosestOrhologSequence(referenceTaxonomy, this.sequences, kegg_taxonomy_ids,
						ncbi_taxonomy_ids, kegg_taxonomy_scores, this.closestOrtholog, orthologsSequences );

				for(int i = 0; i<iterator.size(); i++) {

					String ec_number = iterator.get(i);

					if(!hasLetters(ec_number) && !bypass.contains(ec_number) && ec_number != null) {

						try {

							ConcurrentHashMap<String, AbstractSequence<?>> orthologs = new ConcurrentHashMap<>();

							Map<String, Set<Integer>> genes_ko_modules = new HashMap<>();

							if(gapsIdentification) {

								List<String> kos =	AssembleGPR.getOrthologsByECnumber(ec_number);
								
								for(String ko : kos) {

									List<String> sequenceID = sequenceIdsSet.get(ko);

									if(sequenceID == null || sequenceID.isEmpty()) {

										seq.getOrthologs(ko);

										for(String gene : this.closestOrtholog.get(ko))
											orthologs.put(gene, this.sequences.get(gene));
									}
								}
							}
							else{

								logger.info("Retrieving GPR for {} ...",ec_number);
								AssembleGPR gpr = new AssembleGPR(ec_number);
								Map<String,List<ReactionProteinGeneAssociation>> result = gpr.run();
								logger.info("Retrieved!");

								genes_ko_modules = ModelAPI.loadModule(conn, result);

								logger.info("Genes, KO, modules \t{}",genes_ko_modules);

								for(String ko : genes_ko_modules.keySet()) { 

									List<String> sequenceID = sequenceIdsSet.get(ko);

									if(sequenceID == null || sequenceID.isEmpty()) {

										seq.getOrthologs(ko);

										for(String gene : this.closestOrtholog.get(ko))
											orthologs.put(gene, this.sequences.get(gene));
									}
								}
							}

							logger.info("Orthologs to be searched in genome:\t{}\t{}", orthologs.size(), orthologs.keySet());

							ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet = new ConcurrentLinkedQueue<>();					/////// no outro o concurrent linked queue está aqui

							System.out.println("ORTHOLOGS---->"+orthologs);

							if(orthologs.size()>0) {

								System.out.println("GENOME---->"+this.genome.keySet());
								System.out.println("GENOME SIZE----->"+this.genome.size());

								RunSimilaritySearch search = new RunSimilaritySearch(this.genome, this.similarity_threshold, 
										this.method, orthologs, this.cancel, new AtomicInteger(0), new AtomicInteger(0), AlignmentScoreType.ALIGNMENT);

								search.addObserver(this);
								search.setEc_number(ec_number);

								if(!gapsIdentification)								
									search.setModules(genes_ko_modules);	

								search.setClosestOrthologs(MapUtils.revertMapFromSet(this.closestOrtholog));
								search.setReferenceTaxonomyScore(referenceTaxonomy.size());
								search.setKegg_taxonomy_scores(kegg_taxonomy_scores);
								search.setAnnotatedGenes(this.ecNumbers.get(ec_number));
								search.setReferenceTaxonomyThreshold(this.referenceTaxonomyThreshold);
								search.setCompareToFullGenome(this.compareToFullGenome);

								if(gapsIdentification){
//									boolean recursive = false;
									search.setGapsIdentification(true);
									search.setSubjectFastaFilePath(wsTaxonomyTempFolderPath.concat("GapsFillAnnotationsFile.faa"));
									alignmentContainerSet = search.run_OrthologGapsSearch(sequenceIdsSet, alignmentContainerSet);//, recursive);
								}
								else{
									System.out.println("CHECK2");
									search.setSubjectFastaFilePath(wsTaxonomyTempFolderPath.concat("gprsAnnotationsFile.faa"));
									alignmentContainerSet = search.run_OrthologsSearch(sequenceIdsSet, alignmentContainerSet);
								}


								for (AlignmentCapsule capsule : alignmentContainerSet) {

									if(geneIds.get(capsule.getTarget()) != null) 
										HomologyAPI.loadOrthologsInfo(capsule, geneIds, statement);
								}

								if(gapsIdentification) {

									this.findGapsResult.addAll(alignmentContainerSet);
								}
								else {

									Map<Integer, Set<String>> modules = MapUtils.revertMapFromSet(genes_ko_modules);

									for(int module_id : modules.keySet()) {

										if(search.getSequencesWithoutSimilarities().containsAll(modules.get(module_id)))
											ModelAPI.updateECNumberModule(conn, ec_number, module_id);
										else
											ModelAPI.updateECNumberModule(conn, ec_number, module_id);
									}
								}

							}
							else { 
								
								Map<Integer, Set<String>> modules = MapUtils.revertMapFromSet(genes_ko_modules);

								if(modules.keySet().size()>0) {
									for(int module_id : modules.keySet()) {

										ModelAPI.updateECNumberModule(conn, ec_number, module_id);
									}
								}
							}
						}
						catch (Exception e) {

							ModelAPI.updateECNumberModuleStatus(conn, ec_number, DatabaseProgressStatus.PROCESSING.toString());
							ret = false;
							logger.error("error {}",e.getStackTrace().toString());
							e.printStackTrace();
						}
					}

					if(ret)
						IdentifyGenomeSubunits.setECNumberModuleProcessed(conn, ec_number);

					if(cancel.get())
						i = iterator.size();

					if(this.progress!=null)
						progress.setTime(GregorianCalendar.getInstance().getTimeInMillis() - this.startTime, i, iterator.size());
				}
			}
			conn.closeConnection();
		} 
		catch (Exception e) {

			ret = false;
			logger.error("{}\n{}", e.getMessage(), e);
			e.printStackTrace();
			throw e;
		}
		return ret;
	}



	/**
	 * @return
	 */
	public ConcurrentLinkedQueue<AlignmentCapsule> findGapsResult(){

		return this.findGapsResult;
	}

	/**
	 * Run gene-protein reactions assignment.
	 * 
	 * @param threshold 
	 * @throws SQLException
	 */
	public static Map<String, ReactionsGPR_CI> runGPRsAssignment(double threshold, Connection conn) throws SQLException {


		return ModelAPI.runGPRsAssignment(threshold, conn);
	}

	/**
	 * @param s
	 * @return
	 */
	public static boolean hasLetters(String s) {

		if (s == null) return false;
		for (int i=0; i<s.length(); i++) {
			char c = s.charAt(i);
			if (Character.isLetter(c)) return true;
		}
		return false;
	}



	/**
	 * @return
	 * @throws Exception 
	 */
	public static Map<String, String> getKeggTaxonomyIDs() throws Exception {

		Map<String, String> kegg_taxonomy_ids = new HashMap<>();
		List<String[]> organisms = KeggAPI.getGenomes();

		for(String[] org : organisms)
			kegg_taxonomy_ids.put(org[0], org[1]);

		return kegg_taxonomy_ids;
	}



	/**
	 * @param conn
	 * @param ec_number
	 * @throws SQLException 
	 */
	public static void setECNumberModuleProcessed(Connection conn, String ec_number) throws SQLException {

		ModelAPI.updateECNumberModuleStatus(conn, ec_number,DatabaseProgressStatus.PROCESSED.toString());
	}

	/**
	 * @return the wsTaxonomyTempFolderPath
	 */
	public String getWsTaxonomyTempFolderPath() {
		return wsTaxonomyTempFolderPath;
	}


	/**
	 * @param wsTaxonomyTempFolderPath the wsTaxonomyTempFolderPath to set
	 */
	public void setWsTaxonomyTempFolderPath(String wsTaxonomyTempFolderPath) {
		this.wsTaxonomyTempFolderPath = wsTaxonomyTempFolderPath;
	}


	/**
	 * @param progress
	 */
	public void setProgress(TimeLeftProgress progress) {

		this.progress = progress;
	}

	/**
	 * @param cancel
	 */
	public void setCancel(AtomicBoolean cancel) {

		this.cancel = cancel;
	}


	@Override
	public void update(Observable arg0, Object arg1) {

		setChanged();
		notifyObservers();
	}
}
