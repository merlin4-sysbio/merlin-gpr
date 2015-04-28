package pt.uminho.sysbio.merlin.gpr.rules.core;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.biojava3.core.sequence.ProteinSequence;

import pt.uminho.sysbio.common.bioapis.externalAPI.kegg.KeggAPI;

public class GetClosestOrhologSequence {

	private ConcurrentHashMap<String, ProteinSequence> sequences;
	private ConcurrentHashMap<String, Set<String>> closestOrtholog;
	private Map<String, String> kegg_taxonomy_ids;
	private ConcurrentHashMap<String, Integer> ncbi_taxonomy_ids;
	private List<String> referenceTaxonomy;
	private String ko;
	private ConcurrentHashMap<String, Integer> kegg_taxonomy_scores;
	private ConcurrentHashMap<String, Map<String, List<String>>> orthologsSequences;

	/**
	 * @param ko
	 * @param referenceTaxonomy
	 * @param sequences
	 * @param kegg_taxonomy_ids
	 * @param ncbi_taxonomy_ids
	 * @param kegg_taxonomy_scores
	 * @param closestOrtholog
	 * @param orthologsSequences
	 * @throws Exception
	 */
	public GetClosestOrhologSequence(String ko, List<String> referenceTaxonomy, ConcurrentHashMap<String, ProteinSequence> sequences,
			Map<String, String> kegg_taxonomy_ids, ConcurrentHashMap<String, Integer> ncbi_taxonomy_ids,
			ConcurrentHashMap<String, Integer> kegg_taxonomy_scores, ConcurrentHashMap<String, Set<String>> closestOrtholog, 
			ConcurrentHashMap<String,Map<String,List<String>>> orthologsSequences) throws Exception {

		this.setKo(ko);
		this.setSequences(sequences);
		this.setKegg_taxonomy_ids(kegg_taxonomy_ids);
		this.setNcbi_taxonomy_ids(ncbi_taxonomy_ids);
		this.setClosestOrtholog(closestOrtholog);
		this.referenceTaxonomy = referenceTaxonomy;
		this.setKegg_taxonomy_scores(kegg_taxonomy_scores);
		this.setOrthologsSequences(orthologsSequences);
	}

	/**
	 * @throws Exception
	 */
	public void run() throws Exception {

		this.getOrtholog();
	}

	/**
	 * 
	 * 
	 * @return
	 * @throws Exception
	 */
	public Set<String> getOrtholog() throws Exception {

		if(!this.closestOrtholog.containsKey(ko)) {

			Set<String> genes_ids = this.getGenesForOrtholog();
			this.closestOrtholog.put(ko,genes_ids);

			for(String gene_id : genes_ids) {

				Map<String, List<String>> geneSequences = this.orthologsSequences.get(gene_id);

				String sequence = "";

				for(int i = 1 ; i < geneSequences.get("AASEQ").size(); i++) {

					String aaseq = geneSequences.get("AASEQ").get(i);
					sequence = sequence.concat(aaseq);
				}

				this.sequences.put(gene_id, new ProteinSequence(sequence.toUpperCase()));
			}
		}

		return this.closestOrtholog.get(ko);
	}



	/**
	 * @return 
	 * @throws Exception
	 */
	private Set<String> getGenesForOrtholog() throws Exception {

		String[] findGenes = KeggAPI.findGenes(ko.trim());

		ConcurrentLinkedQueue<String> genes = new ConcurrentLinkedQueue<>();

		for(String gene : findGenes) {

			String[] orgGenes = gene.split(":");
			genes.add(orgGenes[0]);
		}

		SelectClosestOrtholog sco = new SelectClosestOrtholog(this.referenceTaxonomy, this.ncbi_taxonomy_ids, this.kegg_taxonomy_ids, this.kegg_taxonomy_scores, genes);
		sco.run();

		Set<String> gene_id = new HashSet<>();
		int maxScore = -1;
		String maxOrg = "";

		for(String gene : findGenes) {

			String[] orgGenes = gene.split(":");

			String geneOrg = this.kegg_taxonomy_ids.get(orgGenes[0]);

			if(geneOrg!= null && this.ncbi_taxonomy_ids.containsKey(geneOrg)) {

				int score = this.ncbi_taxonomy_ids.get(geneOrg);

				if(this.getMaxScore(score, maxScore)) {

					Map<String, List<String>> geneSequence = KeggAPI.getGenesByID(gene);

					if(geneSequence!=null) {

						gene_id = new HashSet<>();
						maxScore = score;
						maxOrg = orgGenes[0];
						gene_id.add(gene);
						this.orthologsSequences.put(gene, geneSequence);
					}
				}
				else if(maxOrg.equalsIgnoreCase(orgGenes[0])) {

					Map<String, List<String>> geneSequence = KeggAPI.getGenesByID(gene);
					maxScore = score;
					gene_id.add(gene);
					this.orthologsSequences.put(gene, geneSequence);
				}
			}
		}

		return gene_id;
	}


	/**
	 * @param score
	 * @param maxScore
	 * @return
	 */
	private boolean getMaxScore(int score, int maxScore) {

		return score > maxScore;
	}


	/**
	 * @return the ko
	 */
	public String getKo() {
		return ko;
	}


	/**
	 * @param ko the ko to set
	 */
	public void setKo(String ko) {
		this.ko = ko;
	}


	/**
	 * @return the sequences
	 */
	public ConcurrentHashMap<String, ProteinSequence> getSequences() {
		return sequences;
	}


	/**
	 * @param sequences the sequences to set
	 */
	public void setSequences(ConcurrentHashMap<String, ProteinSequence> sequences) {
		this.sequences = sequences;
	}


	/**
	 * @return the ncbi_taxonomy_ids
	 */
	public ConcurrentHashMap<String, Integer> getNcbi_taxonomy_ids() {
		return ncbi_taxonomy_ids;
	}


	/**
	 * @param ncbi_taxonomy_ids the ncbi_taxonomy_ids to set
	 */
	public void setNcbi_taxonomy_ids(ConcurrentHashMap<String, Integer> ncbi_taxonomy_ids) {
		this.ncbi_taxonomy_ids = ncbi_taxonomy_ids;
	}


	/**
	 * @return the kegg_taxonomy_ids
	 */
	public Map<String, String> getKegg_taxonomy_ids() {
		return kegg_taxonomy_ids;
	}


	/**
	 * @param kegg_taxonomy_ids the kegg_taxonomy_ids to set
	 */
	public void setKegg_taxonomy_ids(Map<String, String> kegg_taxonomy_ids) {
		this.kegg_taxonomy_ids = kegg_taxonomy_ids;
	}

	/**
	 * @return the closestOrtholog
	 */
	public ConcurrentHashMap<String, Set<String>> getClosestOrtholog() {
		return closestOrtholog;
	}

	/**
	 * @param closestOrtholog the closestOrtholog to set
	 */
	public void setClosestOrtholog(ConcurrentHashMap<String, Set<String>> closestOrtholog) {
		this.closestOrtholog = closestOrtholog;
	}

	public ConcurrentHashMap<String, Integer> getKegg_taxonomy_scores() {
		return kegg_taxonomy_scores;
	}

	public void setKegg_taxonomy_scores(ConcurrentHashMap<String, Integer> kegg_taxonomy_scores) {
		this.kegg_taxonomy_scores = kegg_taxonomy_scores;
	}

	public ConcurrentHashMap<String, Map<String, List<String>>> getOrthologsSequences() {
		return orthologsSequences;
	}

	public void setOrthologsSequences(ConcurrentHashMap<String, Map<String, List<String>>> orthologsSequences) {
		this.orthologsSequences = orthologsSequences;
	}
}
