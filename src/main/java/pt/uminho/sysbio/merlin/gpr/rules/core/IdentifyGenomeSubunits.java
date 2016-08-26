/**
 * 
 */
package pt.uminho.sysbio.merlin.gpr.rules.core;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Observable;
import java.util.Observer;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.ceb.biosystems.mew.utilities.datastructures.map.MapUtils;
import pt.uminho.sysbio.common.bioapis.externalAPI.kegg.KeggAPI;
import pt.uminho.sysbio.common.bioapis.externalAPI.ncbi.NcbiAPI;
import pt.uminho.sysbio.common.database.connector.databaseAPI.ModelAPI;
import pt.uminho.sysbio.common.database.connector.databaseAPI.containers.gpr.ReactionProteinGeneAssociation;
import pt.uminho.sysbio.common.database.connector.databaseAPI.containers.gpr.ReactionsGPR_CI;
import pt.uminho.sysbio.common.database.connector.datatypes.Connection;
import pt.uminho.sysbio.common.database.connector.datatypes.DatabaseAccess;
import pt.uminho.sysbio.common.local.alignments.core.Enumerators.Method;
import pt.uminho.sysbio.common.local.alignments.core.Enumerators.ThresholdType;
import pt.uminho.sysbio.common.local.alignments.core.PairwiseSequenceAlignement;
import pt.uminho.sysbio.common.local.alignments.core.RunSimilaritySearch;
import pt.uminho.sysbio.merlin.utilities.DatabaseProgressStatus;
import pt.uminho.sysbio.merlin.utilities.TimeLeftProgress;

/**
 * @author ODias
 *
 */
public class IdentifyGenomeSubunits extends Observable implements Observer {

	private static final Logger logger = LoggerFactory.getLogger(IdentifyGenomeSubunits.class);

	private Map<String, List<String>> ecNumbers;
	private Map<String, ProteinSequence> genome;
	private long reference_organism_id;
	private ConcurrentHashMap<String, ProteinSequence> sequences;
	private ConcurrentHashMap<String, Set<String>> closestOrtholog;
	private DatabaseAccess dba;
	private double similarity_threshold;
	private Method method;
	private AtomicBoolean cancel;
	private double referenceTaxonomyThreshold;
	private boolean compareToFullGenome;
	private TimeLeftProgress progress;


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
	public IdentifyGenomeSubunits(Map<String, List<String>> ec_numbers, Map<String, ProteinSequence> genome, long reference_organism_id, 
			DatabaseAccess dba, double similarity_threshold, double referenceTaxonomyThreshold, Method method, 
			boolean compareToFullGenome, AtomicBoolean cancel) {

		this.ecNumbers = ec_numbers;
		this.genome = genome;
		this.reference_organism_id = reference_organism_id;
		this.dba = dba;
		this.similarity_threshold = similarity_threshold;
		this.method = method;
		this.cancel = cancel;
		this.referenceTaxonomyThreshold = referenceTaxonomyThreshold;
		this.compareToFullGenome = compareToFullGenome;
	}


	/**
	 * @throws Exception
	 */
	public boolean runIdentification() throws Exception {

		boolean ret = true;

		try {

			this.sequences = new ConcurrentHashMap<>();
			this.closestOrtholog = new ConcurrentHashMap<>();

			List<String> referenceTaxonomy = NcbiAPI.getReferenceTaxonomy(reference_organism_id);
			logger.info("Reference taxonomy set to {}", referenceTaxonomy);

			ConcurrentHashMap<String, Integer> ncbi_taxonomy_ids = new ConcurrentHashMap<>();
			ConcurrentHashMap<String, Integer> kegg_taxonomy_scores = new ConcurrentHashMap<>();
			ConcurrentHashMap<String, Map<String, List<String>>> orthologsSequences = new ConcurrentHashMap<>();;

			kegg_taxonomy_scores.put("noOrg", 0);
			Map<String, String> kegg_taxonomy_ids = IdentifyGenomeSubunits.getKeggTaxonomyIDs();

			Connection conn = new Connection(this.dba);

			Set<String> bypass = ModelAPI.getECNumbersWithModules(conn);

			long startTime = GregorianCalendar.getInstance().getTimeInMillis();
			List<String> iterator = new ArrayList<>(this.ecNumbers.keySet());

			for(int i = 0; i<iterator.size(); i++) {

				String ec_number = iterator.get(i);

				if(!hasLetters(ec_number) && !bypass.contains(ec_number)  && !this.cancel.get()) {

					try {

						logger.info("Retrieving GPR for {} ...",ec_number);
						AssembleGPR gpr = new AssembleGPR(ec_number);
						Map<String,List<ReactionProteinGeneAssociation>> result = gpr.run();
						logger.info("Retrieved!");

						Map<String, Set<String>> genes_ko_modules = ModelAPI.loadModule(conn, result);
						logger.info("Genes, KO, modules \t{}",genes_ko_modules);

						ConcurrentHashMap<String, ProteinSequence> orthologs = new ConcurrentHashMap<>();
						boolean noOrthologs = true;

						for(String ko : genes_ko_modules.keySet()) { 

							if(!this.cancel.get()){

								List<String> locusTags = ModelAPI.checkDatabase(conn, ko);

								if(locusTags.isEmpty()) {

									GetClosestOrhologSequence seq = new GetClosestOrhologSequence(ko, referenceTaxonomy, this.sequences, kegg_taxonomy_ids,
											ncbi_taxonomy_ids, kegg_taxonomy_scores, this.closestOrtholog, orthologsSequences );
									seq.run();

									for(String gene : this.closestOrtholog.get(ko))
										orthologs.put(gene, this.sequences.get(gene));
								}
								else {

									Map<String, Set<String>> temp = ModelAPI.getOrthologs(ko, conn);

									for(String key : temp.keySet()) {

										noOrthologs = false;

										for(String locus :locusTags) {

											String[] similarityData = new String[4];

											similarityData[0]= key;
											similarityData[1]= locus;
											similarityData[2]= null;
											similarityData[3]= null;

											PairwiseSequenceAlignement.loadOrthologsData(similarityData, conn, ec_number, temp, genes_ko_modules, this.dba.get_database_type());
										}
									}
								}
							}
						}
						logger.info("Orthologs to be searched in genome:\t{}",orthologs.keySet());

						if(orthologs.size()>0 && !this.cancel.get()) {

							RunSimilaritySearch search = new RunSimilaritySearch(this.dba, this.genome, this.similarity_threshold, 
									this.method, orthologs, this.cancel, new AtomicInteger(0), new AtomicInteger(0), -1, ThresholdType.ALIGNMENT);

							//search.addObserver(this);
							search.setEc_number(ec_number);
							search.setModules(genes_ko_modules);
							search.setClosestOrthologs(MapUtils.revertMapFromSet(this.closestOrtholog));
							search.setReferenceTaxonomyScore(referenceTaxonomy.size());
							search.setKegg_taxonomy_scores(kegg_taxonomy_scores);
							search.setAnnotatedGenes(this.ecNumbers.get(ec_number));
							search.setReferenceTaxonomyThreshold(this.referenceTaxonomyThreshold);
							search.setCompareToFullGenome(this.compareToFullGenome);
							search.run_OrthologsSearch();

							Map<String, Set<String>> modules = MapUtils.revertMapFromSet(genes_ko_modules);

							for(String module : modules.keySet()) {

								int module_id = Integer.parseInt(module);

								if(search.getSequencesWithoutSimilarities().containsAll(modules.get(module)))
									ModelAPI.updateECNumberNote(conn, ec_number, module_id, "no_similarities");
							}
						}
						else if(noOrthologs) {

							ModelAPI.updateECNumberNote(conn, ec_number, -1, null);
						}
					}
					catch (Exception e) {

						ModelAPI.updateECNumberStatus(conn, ec_number, DatabaseProgressStatus.PROCESSING.toString());
						ret = false;
						logger.error("error {}",e.getMessage());
					}
				}
				IdentifyGenomeSubunits.setSubunitProcessed(conn, ec_number);

				if(cancel.get())
					i = iterator.size();

				if(this.progress!=null)
					progress.setTime(GregorianCalendar.getInstance().getTimeInMillis() - startTime, i, iterator.size());
			}
			conn.closeConnection();
		} 
		catch (Exception e) {

			ret = false;
			e.printStackTrace();
			logger.error("{}\n{}", e.getMessage(), e);
			throw e;
		}
		return ret;
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
		for(String[] org : organisms) {

			kegg_taxonomy_ids.put(org[0], org[1]);
		}

		return kegg_taxonomy_ids;
	}



	/**
	 * @param conn
	 * @param ec_number
	 * @throws SQLException 
	 */
	public static void setSubunitProcessed(Connection conn, String ec_number) throws SQLException {

		ModelAPI.updateECNumberStatus(conn, ec_number,DatabaseProgressStatus.PROCESSED.toString());
	}

	/**
	 * @param progress
	 */
	public void setProgress(TimeLeftProgress progress) {

		this.progress = progress;
	}


	@Override
	public void update(Observable arg0, Object arg1) {

		setChanged();
		notifyObservers();
	}
}
