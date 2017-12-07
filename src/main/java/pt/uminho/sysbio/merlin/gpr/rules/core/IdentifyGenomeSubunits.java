/**
 * 
 */
package pt.uminho.sysbio.merlin.gpr.rules.core;

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

import pt.uminho.ceb.biosystems.mew.utilities.datastructures.map.MapUtils;
import pt.uminho.sysbio.common.bioapis.externalAPI.kegg.KeggAPI;
import pt.uminho.sysbio.common.bioapis.externalAPI.ncbi.NcbiAPI;
import pt.uminho.sysbio.common.database.connector.databaseAPI.HomologyAPI;
import pt.uminho.sysbio.common.database.connector.databaseAPI.ModelAPI;
import pt.uminho.sysbio.common.database.connector.datatypes.Connection;
import pt.uminho.sysbio.common.database.connector.datatypes.DatabaseAccess;
import pt.uminho.sysbio.common.local.alignments.core.PairwiseSequenceAlignement;
import pt.uminho.sysbio.common.local.alignments.core.RunSimilaritySearch;
import pt.uminho.sysbio.merIin.utilities.capsules.AlignmentCapsule;
import pt.uminho.sysbio.merlin.utilities.DatabaseProgressStatus;
import pt.uminho.sysbio.merlin.utilities.Enumerators.AlignmentScoreType;
import pt.uminho.sysbio.merlin.utilities.Enumerators.Method;
import pt.uminho.sysbio.merlin.utilities.TimeLeftProgress;
import pt.uminho.sysbio.merlin.utilities.containers.gpr.ReactionProteinGeneAssociation;
import pt.uminho.sysbio.merlin.utilities.containers.gpr.ReactionsGPR_CI;

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
	private ConcurrentLinkedQueue<AlignmentCapsule> alignmentResults;



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
	public boolean runIdentification(boolean gapsIdentification) throws Exception {

		boolean ret = true;

		try {

			Connection conn = new Connection(this.dba);
			Statement statement = conn.createStatement();

			this.sequences = new ConcurrentHashMap<>();
			this.closestOrtholog = new ConcurrentHashMap<>();

			List<String> referenceTaxonomy = NcbiAPI.getReferenceTaxonomy(reference_organism_id);
			logger.info("Reference taxonomy set to {}", referenceTaxonomy);

			ConcurrentHashMap<String, Integer> ncbi_taxonomy_ids = new ConcurrentHashMap<>();
			ConcurrentHashMap<String, Integer> kegg_taxonomy_scores = new ConcurrentHashMap<>();
			ConcurrentHashMap<String, Map<String, List<String>>> orthologsSequences = new ConcurrentHashMap<>();;

			kegg_taxonomy_scores.put("noOrg", 0);
			Map<String, String> kegg_taxonomy_ids = IdentifyGenomeSubunits.getKeggTaxonomyIDs();

			Set<String> bypass = null;															//////

			if(!gapsIdentification)														//////// provavelmente pode ser usado pelos gaps
				bypass = ModelAPI.getECNumbersWithModules(conn);						///////// nao e usado nos gaps --> por dentro de condicao

			long startTime = GregorianCalendar.getInstance().getTimeInMillis();

			List<String> iterator = new ArrayList<>(this.ecNumbers.keySet());

			Map<String, Integer> geneIds = ModelAPI.getGeneIds(statement);

			Map<String, List<String>> sequenceIdsSet = ModelAPI.getSequenceIds(statement);

			for(int i = 0; i<iterator.size(); i++) {

				String ec_number = iterator.get(i);

				if(!hasLetters(ec_number) && !bypass.contains(ec_number)  && !this.cancel.get()) {

					try {

						ConcurrentHashMap<String, AbstractSequence<?>> orthologs = new ConcurrentHashMap<>();

						Map<String, Set<String>> genes_ko_modules = new HashMap<>();

						if(!gapsIdentification){

							logger.info("Retrieving GPR for {} ...",ec_number);
							AssembleGPR gpr = new AssembleGPR(ec_number);
							Map<String,List<ReactionProteinGeneAssociation>> result = gpr.run();
							logger.info("Retrieved!");

							genes_ko_modules = ModelAPI.loadModule(conn, result);

							logger.info("Genes, KO, modules \t{}",genes_ko_modules);

							for(String ko : genes_ko_modules.keySet()) { 

								if(!this.cancel.get()){

									List<String> sequenceID = sequenceIdsSet.get(ko);						////////////////

									if(sequenceID == null || sequenceID.isEmpty()) {

										GetClosestOrhologSequence seq = new GetClosestOrhologSequence(ko, referenceTaxonomy, this.sequences, kegg_taxonomy_ids,
												ncbi_taxonomy_ids, kegg_taxonomy_scores, this.closestOrtholog, orthologsSequences );

										seq.run();

										for(String gene : this.closestOrtholog.get(ko))
											orthologs.put(gene, this.sequences.get(gene));
									}
								}
							}
						}
						else{
							List<String> kos =	AssembleGPR.getOrthologsByECnumber(ec_number);

							for(String ko : kos) {

								List<String> sequenceID = sequenceIdsSet.get(ko);

								if(sequenceID == null || sequenceID.isEmpty()) {

									GetClosestOrhologSequence seq = new GetClosestOrhologSequence(ko, referenceTaxonomy, this.sequences, kegg_taxonomy_ids,
											ncbi_taxonomy_ids, kegg_taxonomy_scores, this.closestOrtholog, orthologsSequences );
									seq.run();

									Set<String> orthologsID = this.closestOrtholog.get(ko);

									for(String gene : orthologsID)
										orthologs.put(gene, this.sequences.get(gene));
								}
							}
						}


						logger.info("Orthologs to be searched in genome:\t{}",orthologs.keySet());

						ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet = new ConcurrentLinkedQueue<>();					/////// no outro o concurrent linked queue está aqui

						if(orthologs.size()>0 && !this.cancel.get()) {

							RunSimilaritySearch search = new RunSimilaritySearch(this.genome, this.similarity_threshold, 
									this.method, orthologs, this.cancel, new AtomicInteger(0), new AtomicInteger(0), AlignmentScoreType.ALIGNMENT);

							//search.addObserver(this);
							search.setEc_number(ec_number);
							if(!gapsIdentification)
								search.setModules(genes_ko_modules);											
							search.setClosestOrthologs(MapUtils.revertMapFromSet(this.closestOrtholog));
							search.setReferenceTaxonomyScore(referenceTaxonomy.size());
							search.setKegg_taxonomy_scores(kegg_taxonomy_scores);
							search.setAnnotatedGenes(this.ecNumbers.get(ec_number));
							search.setReferenceTaxonomyThreshold(this.referenceTaxonomyThreshold);
							search.setCompareToFullGenome(this.compareToFullGenome);

							if(!gapsIdentification)
								alignmentContainerSet = search.run_OrthologsSearch(sequenceIdsSet, alignmentContainerSet);			/////// aqui e usado outro método!!! condicao
							else
								alignmentContainerSet = search.run_OrthologGapsSearch(sequenceIdsSet, alignmentContainerSet);

							for (AlignmentCapsule capsule : alignmentContainerSet) {

								if(geneIds.get(capsule.getTarget()) != null) 
									HomologyAPI.loadOrthologsInfo(capsule, geneIds, statement);
							}

							if(!gapsIdentification){
								Map<String, Set<String>> modules = MapUtils.revertMapFromSet(genes_ko_modules);					///////// esta parte nao e usada no outro

								for(String module : modules.keySet()) {

									int module_id = Integer.parseInt(module);

									if(search.getSequencesWithoutSimilarities().containsAll(modules.get(module)))
										ModelAPI.updateECNumberNote(conn, ec_number, module_id, "no_similarities");
								}
							}
						}
						else if(orthologs.size() == 0) {

							ModelAPI.updateECNumberNote(conn, ec_number, -1, null);								////// aqui -1??
						}

					}
					catch (Exception e) {

						ModelAPI.updateECNumberStatus(conn, ec_number, DatabaseProgressStatus.PROCESSING.toString());
						ret = false;
						logger.error("error {}",e.getStackTrace().toString());
						e.printStackTrace();
					}
				}
				if(ret)
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
			logger.error("{}\n{}", e.getMessage(), e);
			e.printStackTrace();
			throw e;
		}
		return ret;
	}
	
	
	


//	/**
//	 * @throws Exception
//	 */
//	public boolean runIdentification() throws Exception {
//
//		boolean ret = true;
//
//		try {
//			
//			Connection conn = new Connection(this.dba);
//			Statement statement = conn.createStatement();
//
//			this.sequences = new ConcurrentHashMap<>();
//			this.closestOrtholog = new ConcurrentHashMap<>();
//
//			List<String> referenceTaxonomy = NcbiAPI.getReferenceTaxonomy(reference_organism_id);
//			logger.info("Reference taxonomy set to {}", referenceTaxonomy);
//
//			ConcurrentHashMap<String, Integer> ncbi_taxonomy_ids = new ConcurrentHashMap<>();
//			ConcurrentHashMap<String, Integer> kegg_taxonomy_scores = new ConcurrentHashMap<>();
//			ConcurrentHashMap<String, Map<String, List<String>>> orthologsSequences = new ConcurrentHashMap<>();;
//
//			kegg_taxonomy_scores.put("noOrg", 0);
//			Map<String, String> kegg_taxonomy_ids = IdentifyGenomeSubunits.getKeggTaxonomyIDs();
//
//			Set<String> bypass = ModelAPI.getECNumbersWithModules(conn);						///////// nao e usado nos gaps --> por dentro de condicao
//
//			long startTime = GregorianCalendar.getInstance().getTimeInMillis();
//			
//			List<String> iterator = new ArrayList<>(this.ecNumbers.keySet());
//			
//			Map<String, Integer> geneIds = ModelAPI.getGeneIds(statement);
//			
//			Map<String, List<String>> sequenceIdsSet = ModelAPI.getSequenceIds(statement);
//			
//			for(int i = 0; i<iterator.size(); i++) {
//																								//////// falta aqui uma linha dos gaps ---> condicao
//				String ec_number = iterator.get(i);
//
//				if(!hasLetters(ec_number) && !bypass.contains(ec_number)  && !this.cancel.get()) {
//
//					try {
//																								/////////////////////////////////////
//						logger.info("Retrieving GPR for {} ...",ec_number);
//						AssembleGPR gpr = new AssembleGPR(ec_number);
//						Map<String,List<ReactionProteinGeneAssociation>> result = gpr.run();
//						logger.info("Retrieved!");
//
//						Map<String, Set<String>> genes_ko_modules = ModelAPI.loadModule(conn, result);
//						
//						logger.info("Genes, KO, modules \t{}",genes_ko_modules);
//
//						ConcurrentHashMap<String, AbstractSequence<?>> orthologs = new ConcurrentHashMap<>();
//						
//						for(String ko : genes_ko_modules.keySet()) { 
//
//							if(!this.cancel.get()){
//								
//								List<String> sequenceID = sequenceIdsSet.get(ko);						////////////////
//								
//								if(sequenceID == null || sequenceID.isEmpty()) {
//
//									GetClosestOrhologSequence seq = new GetClosestOrhologSequence(ko, referenceTaxonomy, this.sequences, kegg_taxonomy_ids,
//											ncbi_taxonomy_ids, kegg_taxonomy_scores, this.closestOrtholog, orthologsSequences );
//
//									seq.run();
//
//									for(String gene : this.closestOrtholog.get(ko))
//										orthologs.put(gene, this.sequences.get(gene));
//								}
//							}
//						}
//						logger.info("Orthologs to be searched in genome:\t{}",orthologs.keySet());
//																											/////// no outro o concurrent linked queue esta aqui
//						if(orthologs.size()>0 && !this.cancel.get()) {
//
//							RunSimilaritySearch search = new RunSimilaritySearch(this.genome, this.similarity_threshold, 
//									this.method, orthologs, this.cancel, new AtomicInteger(0), new AtomicInteger(0), AlignmentScoreType.ALIGNMENT);
//
//							//search.addObserver(this);
//							search.setEc_number(ec_number);
//							search.setModules(genes_ko_modules);											///// o outro nao passa mas acho que nao ha problema em passar
//							search.setClosestOrthologs(MapUtils.revertMapFromSet(this.closestOrtholog));
//							search.setReferenceTaxonomyScore(referenceTaxonomy.size());
//							search.setKegg_taxonomy_scores(kegg_taxonomy_scores);
//							search.setAnnotatedGenes(this.ecNumbers.get(ec_number));
//							search.setReferenceTaxonomyThreshold(this.referenceTaxonomyThreshold);
//							search.setCompareToFullGenome(this.compareToFullGenome);
//							
//							ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet = new ConcurrentLinkedQueue<>();
//							
//							alignmentContainerSet = search.run_OrthologsSearch(sequenceIdsSet, alignmentContainerSet);			/////// aqui e usado outro metodo!!! condicao
//							
//							for (AlignmentCapsule capsule : alignmentContainerSet) {
//								
//								if(geneIds.get(capsule.getTarget()) != null) 
//									HomologyAPI.loadOrthologsInfo(capsule, geneIds, statement);
//							}
//							
//							Map<String, Set<String>> modules = MapUtils.revertMapFromSet(genes_ko_modules);					///////// esta parte nao e usada no outro
//							
//							for(String module : modules.keySet()) {
//								
//								int module_id = Integer.parseInt(module);
//								
//								if(search.getSequencesWithoutSimilarities().containsAll(modules.get(module)))
//									ModelAPI.updateECNumberNote(conn, ec_number, module_id, "no_similarities");
//							}
//						}
//						else if(orthologs.size() == 0) {
//
//							ModelAPI.updateECNumberNote(conn, ec_number, -1, null);
//						}
//						
//					}
//					catch (Exception e) {
//
//						ModelAPI.updateECNumberStatus(conn, ec_number, DatabaseProgressStatus.PROCESSING.toString());
//						ret = false;
//						logger.error("error {}",e.getStackTrace().toString());
//						e.printStackTrace();
//					}
//				}
//				if(ret)
//					IdentifyGenomeSubunits.setSubunitProcessed(conn, ec_number);
//
//				if(cancel.get())
//					i = iterator.size();
//
//				if(this.progress!=null)
//					progress.setTime(GregorianCalendar.getInstance().getTimeInMillis() - startTime, i, iterator.size());
//			}
//			conn.closeConnection();
//		} 
//		catch (Exception e) {
//
//			ret = false;
//			logger.error("{}\n{}", e.getMessage(), e);
//			e.printStackTrace();
//			throw e;
//		}
//		return ret;
//	}
//
//	/**
//	 * @throws Exception
//	 */
//	public boolean runGapsIdentification() throws Exception {
//
////		this.alignmentResults = alignmentResults;
//		boolean ret = true;
//
//		try {
//			
//			Connection conn = new Connection(this.dba);
//			Statement statement = conn.createStatement();
//
//			this.sequences = new ConcurrentHashMap<>();
//			this.closestOrtholog = new ConcurrentHashMap<>();
//
//			List<String> referenceTaxonomy = NcbiAPI.getReferenceTaxonomy(reference_organism_id);
//			logger.info("Reference taxonomy set to {}", referenceTaxonomy);
//
//			ConcurrentHashMap<String, Integer> ncbi_taxonomy_ids = new ConcurrentHashMap<>();
//			ConcurrentHashMap<String, Integer> kegg_taxonomy_scores = new ConcurrentHashMap<>();
//			ConcurrentHashMap<String, Map<String, List<String>>> orthologsSequences = new ConcurrentHashMap<>();;
//
//			kegg_taxonomy_scores.put("noOrg", 0);
//			Map<String, String> kegg_taxonomy_ids = IdentifyGenomeSubunits.getKeggTaxonomyIDs();
//
//			long startTime = GregorianCalendar.getInstance().getTimeInMillis();
//
//			List<String> iterator = new ArrayList<>(this.ecNumbers.keySet());
//			
//			Map<String, Integer> geneIds = ModelAPI.getGeneIds(statement);
//			
//			Map<String, List<String>> sequenceIdsSet = ModelAPI.getSequenceIds(statement);
//			
//			for(int i = 0; i<iterator.size(); i++) {
//				
//				String ec_number = iterator.get(i);
//
//				List<String> kos =	AssembleGPR.getOrthologsByECnumber(ec_number);							/////////////////// Nao e usado no outro
//
//				if(!this.cancel.get()) {
//
//					try {
//
//						ConcurrentHashMap<String, AbstractSequence<?>> orthologs = new ConcurrentHashMap<>();
//
//						if(!this.cancel.get()){
//
//							
//						}
//						logger.info("Orthologs to be searched in genome:\t{}",orthologs.keySet());
//
//						ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet = new ConcurrentLinkedQueue<>();
//						
//						if(orthologs.size()>0 && !this.cancel.get()) {
//
//							RunSimilaritySearch search = new RunSimilaritySearch(this.genome, this.similarity_threshold, 
//									this.method, orthologs, this.cancel, new AtomicInteger(0), new AtomicInteger(0), AlignmentScoreType.ALIGNMENT);
//
//							//search.addObserver(this);
//							search.setEc_number(ec_number);									//////	tambem pode passar os KOS acho eu
//							search.setClosestOrthologs(MapUtils.revertMapFromSet(this.closestOrtholog));
//							search.setReferenceTaxonomyScore(referenceTaxonomy.size());
//							search.setKegg_taxonomy_scores(kegg_taxonomy_scores);
//							search.setAnnotatedGenes(this.ecNumbers.get(ec_number));
//							search.setReferenceTaxonomyThreshold(this.referenceTaxonomyThreshold);
//							search.setCompareToFullGenome(this.compareToFullGenome);
//							
//							alignmentContainerSet = search.run_OrthologGapsSearch(sequenceIdsSet, alignmentContainerSet);
//						}
//						
//						for (AlignmentCapsule capsule : alignmentContainerSet){
//							
//							if(geneIds.get(capsule.getTarget()) != null) 
//								HomologyAPI.loadOrthologsInfo(capsule, geneIds, statement);
//						}
//
////						for(String ko :  alignmentResults.keySet())
////							System.out.println(ko+"\t"+alignmentResults.get(ko)[0]+"\t"+alignmentResults.get(ko)[1]+"\t"+alignmentResults.get(ko)[2]+"\t"+alignmentResults.get(ko)[3]);
//
//					}
//					catch (Exception e) {
//
//						ModelAPI.updateECNumberStatus(conn, ec_number, DatabaseProgressStatus.PROCESSING.toString());
//						ret = false;
//						logger.error("error {}",e.getMessage());
//					}
//				}
//				if(ret)
//					IdentifyGenomeSubunits.setSubunitProcessed(conn, ec_number);
//				
//				if(cancel.get())
//					i = iterator.size();
//
//				if(this.progress!=null)
//					progress.setTime(GregorianCalendar.getInstance().getTimeInMillis() - startTime, i, iterator.size());
//			}
//			conn.closeConnection();
//		} 
//		catch (Exception e) {
//
//			ret = false;
//			e.printStackTrace();
//			logger.error("{}\n{}", e.getMessage(), e);
//			throw e;
//		}
//		return ret;	
//	}

	public ConcurrentLinkedQueue<AlignmentCapsule> getAlignmentResults(){
		
		return alignmentResults;
			
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
