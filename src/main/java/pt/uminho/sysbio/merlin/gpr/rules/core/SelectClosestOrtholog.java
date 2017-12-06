package pt.uminho.sysbio.merlin.gpr.rules.core;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;

import pt.uminho.ceb.biosystems.mew.utilities.datastructures.map.MapUtils;
import pt.uminho.sysbio.common.bioapis.externalAPI.ncbi.NcbiAPI;

/**
 * @author ODias
 *
 */
public class SelectClosestOrtholog //implements Runnable 
{

	private Map<String, String> kegg_taxonomy_ids;
	private ConcurrentHashMap<String, Integer> ncbi_taxonomy_ids;
	private List<String> referenceTaxonomy;
	private ConcurrentLinkedQueue<String> genes;
	private ConcurrentHashMap<String, Integer> kegg_taxonomy_scores;
	private Map<String, Set<String>> kegg_taxonomy_reverse;

	/**
	 * @param referenceTaxonomy
	 * @param ncbi_taxonomy_ids
	 * @param kegg_taxonomy_ids
	 * @param kegg_taxonomy_scores 
	 * @param genes
	 * @throws Exception
	 */
	public SelectClosestOrtholog(List<String> referenceTaxonomy, ConcurrentHashMap<String, Integer> ncbi_taxonomy_ids, Map<String, String> kegg_taxonomy_ids,
			ConcurrentHashMap<String, Integer> kegg_taxonomy_scores, ConcurrentLinkedQueue<String> genes) throws Exception {

		this.ncbi_taxonomy_ids = ncbi_taxonomy_ids;
		this.referenceTaxonomy = referenceTaxonomy;
		this.kegg_taxonomy_ids = kegg_taxonomy_ids;
		this.genes = genes;
		this.kegg_taxonomy_scores = kegg_taxonomy_scores;
		this.kegg_taxonomy_reverse = MapUtils.revertMap(this.kegg_taxonomy_ids);
		//this.processTaxId(referenceTaxonomy, null);
	}

	/**
	 * @throws Exception
	 */
	public void run() throws Exception {

		List<String> genesOrg  = new ArrayList<>();

//		for(String gene : this.genes)
//			genesOrg.add(this.kegg_taxonomy_ids.get(gene.split(":")[0]));

		for(String gene : this.genes){
			
			String taxID = this.kegg_taxonomy_ids.get(gene.split(":")[0]);
			
			if (taxID != null && !taxID.equals(""))
				genesOrg.add(taxID);
		}
		
		this.processTaxIds(genesOrg);
	}

	/**
	 * @param tax_ids
	 * @throws Exception
	 */
	private void processTaxIds(List<String> tax_ids_list) throws Exception {
		
		List<List<String>> lists = this.splitLists(tax_ids_list, 100);

		for(List<String> tax_ids : lists) {

			Map<String,String[]> ncbi_ids = NcbiAPI.getTaxonList(tax_ids.toString());

			for(String tax_id : tax_ids) {

				int score = -1;

				if(ncbi_ids.containsKey(tax_id)) {

					String[] taxonomy = ncbi_ids.get(tax_id)[1].split(";");

					List<String> match = new ArrayList<>();
					
					for(String t : taxonomy)
						match.add(t.trim());

					match.add(ncbi_ids.get(tax_id)[0].trim());
					match.retainAll(this.referenceTaxonomy);

					score = match.size();

					this.ncbi_taxonomy_ids.put(tax_id, score);

					for(String org : this.kegg_taxonomy_reverse.get(tax_id))					
						this.kegg_taxonomy_scores.put(org, score);
				}
			}
		}
	}

	/**
	 * @param list
	 * @param dimension
	 * @return
	 */
	public List<List<String>> splitLists(List<String> list, int dimension) {

		List<List<String>> lists = new ArrayList<>();

		if(list.size()>dimension) {

			lists.add(list.subList(0, dimension));
			lists.addAll(this.splitLists(list.subList(dimension, list.size()), dimension));
		}
		else {

			lists.add(list);
		}

		return lists;
	}

	//	@Override
	//	public void run() {
	//
	//		while(!this.genes.isEmpty()) {
	//
	//			String gene = this.genes.poll();
	//
	//			try {
	//
	//				this.getGenesForOrtholog(gene);
	//			}
	//			catch (Exception e) {
	//
	//				e.printStackTrace();
	//				this.genes.offer(gene);
	//			}
	//		}
	//	}
	//	
	//	
	//	/**
	//	 * @param gene
	//	 * @return
	//	 * @throws Exception
	//	 */
	//	private void getGenesForOrtholog(String gene) throws Exception {
	//
	//		String[] orgGenes = gene.split(":");
	//
	//		if(!this.ncbi_taxonomy_ids.keySet().contains(orgGenes[0])) {
	//
	//			String[] tax_data = KeggAPI.findGenome(this.kegg_taxonomy_ids.get(orgGenes[0]));
	//
	//			Pattern pattern = Pattern.compile("(\\d{3,7})");
	//			Matcher matcher = pattern.matcher(tax_data[0]);
	//
	//			if (matcher.find()) {
	//
	//				String tax_id = matcher.group();
	//
	//				this.processTaxId(tax_id, orgGenes[0]);
	//			}
	//		}
	//	}
	//	/**
	//	 * 
	//	 * 
	//	 * @param tax_id
	//	 * @param tax_code
	//	 * @throws Exception
	//	 */
	//	private void processTaxId(String tax_id, String tax_code) throws Exception {
	//
	//		int score = -1;
	//		if(tax_code==null) {
	//
	//			if(this.referenceTaxonomy.isEmpty()) {
	//
	//				NCBI_Taxon_Stub_API ncsa = new NCBI_Taxon_Stub_API(10);
	//				Map<String,String[]> ncbi_ids= ncsa.getTaxonList(tax_id, 0);
	//
	//				String[] taxonomy = ncbi_ids.get(tax_id)[1].split(";");
	//
	//				this.referenceTaxonomy = new ConcurrentLinkedQueue<>();
	//
	//				for(String t : taxonomy) {
	//
	//					referenceTaxonomy.add(t.trim());
	//				}
	//
	//				this.referenceTaxonomy.add(ncbi_ids.get(tax_id)[0].trim());
	//
	//				score = this.referenceTaxonomy.size();
	//				this.ncbi_taxonomy_ids.put(tax_id, score);
	//			}
	//
	//			//System.out.println(referenceTaxonomy);
	//			//System.out.println(tax_id+"\t"+tax_code+"\t"+this.ncbi_taxonomy_ids.get(tax_id));
	//
	//		}
	//		else {
	//
	//			NCBI_Taxon_Stub_API ncsa = new NCBI_Taxon_Stub_API(10);
	//			Map<String,String[]> ncbi_ids= ncsa.getTaxonList(tax_id, 0);
	//
	//			String[] taxonomy = ncbi_ids.get(tax_id)[1].split(";");
	//			List<String> match = new ArrayList<>();
	//
	//			for(String t : taxonomy) {
	//
	//				match.add(t.trim());
	//			}
	//
	//			match.add(ncbi_ids.get(tax_id)[0].trim());
	//			match.retainAll(this.referenceTaxonomy);
	//
	//			score = match.size();
	//			this.ncbi_taxonomy_ids.put(tax_code, score);
	//
	//			//System.out.println(tax_id+"\t"+tax_code+"\t"+this.ncbi_taxonomy_ids.get(tax_code));
	//		}
	//	}

}
