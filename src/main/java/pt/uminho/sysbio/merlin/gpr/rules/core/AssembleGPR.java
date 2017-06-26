package pt.uminho.sysbio.merlin.gpr.rules.core;

import java.io.ByteArrayInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.sysbio.common.bioapis.externalAPI.kegg.KeggAPI;
import pt.uminho.sysbio.common.database.connector.databaseAPI.containers.gpr.GeneAssociation;
import pt.uminho.sysbio.common.database.connector.databaseAPI.containers.gpr.ModuleCI;
import pt.uminho.sysbio.common.database.connector.databaseAPI.containers.gpr.ProteinGeneAssociation;
import pt.uminho.sysbio.common.database.connector.databaseAPI.containers.gpr.ReactionProteinGeneAssociation;
import pt.uminho.sysbio.common.database.connector.datatypes.Enumerators.ModuleType;
import pt.uminho.sysbio.merlin.gpr.rules.core.input.KeggModulesParser;
import pt.uminho.sysbio.merlin.gpr.rules.grammar.KEGGOrthologyParser;


/**
 * @author ODias
 *
 */
public class AssembleGPR {

	private static final Logger logger = LoggerFactory.getLogger(AssembleGPR.class);

	private String ec_number;
	private Map<String, List<String>> orthologsReactions;
	private List<String> orthologsEnzymes;
	private List<String> reactionsEnzymes;
	private boolean isPartialECnumber;

	//TODO vero que fazer com o menos!!! \\+ \\-

	/**
	 * @param ec_number
	 * @throws IOException 
	 * @throws  
	 */
	public AssembleGPR (String ec_number) throws Exception {

		this.ec_number = ec_number;
		this.isPartialECnumber = ec_number.contains(".-");
		this.orthologsReactions = new HashMap<String, List<String>>();
		//this.initialOrthologsReactions = new HashMap<String, List<String>>();
		this.orthologsEnzymes = new ArrayList<String>();
		this.reactionsEnzymes  = new ArrayList<String>();

		//System.out.println("search ortholog by ec: "+"\n"+"rest.kegg.jp/link/ko/"+ec_number);
		//System.out.println("search reaction by ec: "+"\n"+"rest.kegg.jp/link/reaction/"+ec_number);
	}

	/**
	 * @return
	 * @throws Exception
	 */
	public Map<String,List<ReactionProteinGeneAssociation>> run() throws Exception {

		this.orthologsEnzymes.addAll(AssembleGPR.getOrthologsByECnumber(this.ec_number));

		if(this.isPartialECnumber)
			this.getReactionsByOrthologue();
		else
			this.getReactionsByECnumber();

		this.filterOrthologyByReactions();

		Map<String,List<ReactionProteinGeneAssociation>> result = this.getModules();

		return result;
	}

	/**
	 * @return
	 * @throws Exception
	 */
	public static List<String> getOrthologsByECnumber(String ec_number) throws Exception {

		List<String> orthologs_temp = KeggAPI.getOrthologsByECnumber(ec_number);
		List<String> orthologsEnzymes = new ArrayList<>();
		
		for(String ko : orthologs_temp) {

			StringTokenizer sTok = new StringTokenizer(ko,"ko:");

			while (sTok.hasMoreTokens()){
				String tok = sTok.nextToken();
				Pattern pattern = Pattern.compile("(K\\d{5})");
				Matcher matcher = pattern.matcher(tok);
				if (matcher.find()) {

					orthologsEnzymes.add(matcher.group());
				}
			}
		}

		return orthologsEnzymes;
	}

	/**
	 * @return
	 * @throws Exception
	 */
	public List<String> getReactionsByOrthologue() throws Exception {

		for(String orthologue : this.orthologsEnzymes) {

			Set<String> reactions_temp = KeggAPI.getReactionsByOrthology(orthologue);

			for(String r : reactions_temp) {

				Pattern pattern2 = Pattern.compile("(R\\d{5})");

				Matcher matcher2 = pattern2.matcher(r);
				if (matcher2.find() && !this.reactionsEnzymes.contains(matcher2.group())) {

					this.reactionsEnzymes.add(matcher2.group());
				}
			}

		}
		return this.reactionsEnzymes;
	}

	/**
	 * @return
	 * @throws Exception
	 */
	public List<String> getReactionsByECnumber () throws Exception {

		Set<String> reactions_temp = KeggAPI.getReactionsByEnzymes(ec_number);

		for(String r : reactions_temp) {

			Pattern pattern2 = Pattern.compile("(R\\d{5})");

			Matcher matcher2 = pattern2.matcher(r);
			if (matcher2.find()) {

				this.reactionsEnzymes.add(matcher2.group());
			}
		}
		return this.reactionsEnzymes;
	}

	/**
	 * @return
	 * @throws Exception
	 */
	public Map<String, List<String>> filterOrthologyByReactions() throws Exception {

		for(int j=0;j<this.reactionsEnzymes.size();j++) {

			String reaction = this.reactionsEnzymes.get(j);

			List<String> reactionOrthologs = this.getOrthologsByReaction(reaction);

			List<String> reactionOrthologsClone = new ArrayList<>();
			reactionOrthologsClone.addAll(reactionOrthologs);
			reactionOrthologsClone.retainAll(this.orthologsEnzymes);

			if(reactionOrthologsClone.size()>0) {

				//this.initialOrthologsReactions.put(reaction, new ArrayList<>(reactionOrthologs));
				reactionOrthologs.removeAll(this.orthologsEnzymes);

				if(reactionOrthologs.size()>0)
					reactionOrthologs = this.getOrthologsWihtoutECnumber(reactionOrthologs);

				reactionOrthologs.addAll(reactionOrthologsClone);

				//System.out.println("Reaction "+reaction+" List of orthologs of interest:\t"+reactionOrthologs);

				this.orthologsReactions.put(reaction, reactionOrthologs);
			}
		}

		return orthologsReactions;
	}


	/**
	 * @param reaction
	 * @return
	 * @throws Exception
	 */
	private List<String> getOrthologsByReaction(String reaction) throws Exception {

		List<String> orthologs_temp = KeggAPI.getOrthologsByReaction(reaction);

		List<String> orthologs = new ArrayList<String>();

		for(String ortholog : orthologs_temp) {

			StringTokenizer sTokm = new StringTokenizer(ortholog,"ko:");  

			while (sTokm.hasMoreTokens()) {

				String tokm = sTokm.nextToken();
				Pattern patternm = Pattern.compile("(K\\d{5})");
				Matcher matcherm = patternm.matcher(tokm);
				if (matcherm.find()) {

					if(!orthologs.contains(matcherm.group()))
						orthologs.add(matcherm.group());
				}
			}
		}
		return orthologs;
	}

	/**
	 * @return
	 * @throws Exception
	 */
	public Map<String, List<ReactionProteinGeneAssociation>> getModules() throws Exception {

		Map<String, List<ReactionProteinGeneAssociation>> gpr_rules = new HashMap<>();

		for(String reaction : this.orthologsReactions.keySet()) {

			String procura = this.ec_number+"+"+reaction;

			for (String ortholog : this.orthologsReactions.get(reaction)) {

				procura +="+"+ortholog;
			}

			//System.out.println("Cross-reference module by ec, reaction and orthologs:\t rest.kegg.jp/link/module/"+procura);
			String modulesQuery = KeggAPI.getModulesStringByQuery(procura);

			List<String> modules = this.parseModules(modulesQuery);

			List<ReactionProteinGeneAssociation> rpg = this.verifyModules(modules, reaction, modulesQuery);

			if(rpg != null && rpg.size()>0)
				gpr_rules.put(reaction, rpg);

		}

		return gpr_rules;
	}


	/**
	 * @param modules
	 * @return
	 * @throws Exception
	 */
	private List<String> parseModules(String modules) throws Exception{

		String[] rows = modules.split("\n");
		List<String> returnModules = new ArrayList<String>();

		for(String row : rows) {

			StringTokenizer sTokmod = new StringTokenizer(row,"md:");  
			while (sTokmod.hasMoreTokens()){

				String tokmod = sTokmod.nextToken();
				Pattern patternmod = Pattern.compile("(M\\d{5})");
				Matcher matchermod = patternmod.matcher(tokmod);

				if (matchermod.find()) {

					if(!returnModules.contains(matchermod.group())) {

						returnModules.add(matchermod.group());
					}
				}
			}
		}
		return returnModules;
	}


	/**
	 * @param modules
	 * @param reaction
	 * @param modulesQuery
	 * @return
	 * @throws Exception
	 */
	private List<ReactionProteinGeneAssociation> verifyModules(List<String> modules, String reaction, String modulesQuery) throws Exception{

		List<ReactionProteinGeneAssociation> gpr_list = new ArrayList<>();

		ReactionProteinGeneAssociation gpr = this.verifyModule(modules, reaction, modulesQuery);

		if(gpr!=null)
			gpr_list.add(gpr);

		return gpr_list;
	}

	/**
	 * @param module
	 * @param reaction
	 * @param reactionOrthologs
	 * @throws Exception
	 */
	private ReactionProteinGeneAssociation verifyModule(List<String> modules, String reaction, String modulesQuery) throws Exception {

		ReactionProteinGeneAssociation gpr = new ReactionProteinGeneAssociation(reaction);
		ProteinGeneAssociation protein_rule = new ProteinGeneAssociation(this.ec_number);

		for(String module : modules) {

			ModuleConfirmation moduleConfirmation = new ModuleConfirmation(module);

			String[] rows = modulesQuery.split("\n");

			for(String row : rows) {

				Pattern patternbmodul = Pattern.compile(module);
				Matcher matcherbmodul = patternbmodul.matcher(row);

				if(matcherbmodul.find()) {

					Pattern patternc = Pattern.compile("ec:");
					Matcher matcherc = patternc.matcher(row);
					Pattern patternd = Pattern.compile("rn:");
					Matcher matcherd = patternd.matcher(row);
					Pattern patternf = Pattern.compile("K\\d{5}");
					Matcher matcherf = patternf.matcher(row);

					boolean ec_match = matcherc.find();
					boolean reaction_match = matcherd.find();
					boolean orthologs_match = matcherf.find();

					if(ec_match){

						moduleConfirmation.setHasECnumber(true);
					}
					else if(reaction_match) {

						moduleConfirmation.setHasReaction(true);
					}
					else if(orthologs_match) {

						moduleConfirmation.addOrtholog(matcherf.group());
					}
				}
			}

			if(this.isPartialECnumber)
				moduleConfirmation.setHasECnumber(true);

			//System.out.println("Module\t"+module+"\t for reaction\t"+reaction+"\t:\t"+moduleConfirmation);
			ModuleType moduleType = this.confirmModule(moduleConfirmation, reaction); 

			if(moduleType != null) {

				String s = KeggAPI.getModuleEntry(module);

				KeggModulesParser k = new KeggModulesParser(s);

				ModuleCI mic = this.parseModuleCI(k, module, moduleType);

				Set<String> orthologsOfInterest = this.getOrthologsOfInterest(k, moduleType, reaction, new ArrayList<String> (moduleConfirmation.getOrthologs())); 

				if(orthologsOfInterest.size()>0) {

					List<GeneAssociation> geneAssociationList = this.getdefinition(orthologsOfInterest, reaction, module, moduleType, mic);

					if(geneAssociationList!=null)
						protein_rule.addAllGeneAssociation(geneAssociationList);
				}
			}
		}

		gpr.addProteinGeneAssociation(protein_rule);

		if(gpr.getProteinGeneAssociation().get(this.ec_number).getGenes().isEmpty())
			return null;

		return gpr;
	}

	/**
	 * 
	 * @param moduleConfirmation
	 * @param reaction
	 * @return
	 */
	private ModuleType confirmModule(ModuleConfirmation moduleConfirmation, String reaction){

		if(!moduleConfirmation.isHasECnumber() || 
				moduleConfirmation.getOrthologs() == null || 
				orthologsReactions.get(reaction) == null)
			return null;

		if(orthologsReactions.get(reaction).containsAll(moduleConfirmation.getOrthologs())) {

			if(moduleConfirmation.isHasReaction())				
				return ModuleType.Pathway;
			else
				return ModuleType.Complex;


		}
		return null;
	}

	/**
	 * 
	 * 
	 * @param moduleEntry
	 * @param moduleType
	 * @param reaction
	 * @param orthologs
	 * @return
	 * @throws Exception
	 */
	private Set<String> getOrthologsOfInterest(KeggModulesParser modules, ModuleType moduleType, String reaction, List<String> orthologs) throws Exception {

		Set<String> ret = new HashSet<>();

		List<String> expregene = new ArrayList<>();
		expregene.addAll(orthologs);
		expregene.add(ec_number);
		if(moduleType.equals(ModuleType.Pathway))
			expregene.add(reaction);
		expregene.add(moduleType.toString());

		if(modules.getEntry().contains(moduleType.toString()))
			expregene.remove(moduleType.toString());

		for(String line : modules.getOrthology().split("\n")) {

			if(line.contains(ec_number))
				expregene.remove(ec_number);

			if(line.contains(reaction))
				expregene.remove(reaction);

			for(String ortholog : orthologs) {

				if(line.contains(ortholog)) {

					expregene.remove(ortholog);
				}
				ret.add(ortholog);
			}

			if(expregene.isEmpty())
				return ret;

			if(moduleType.equals(ModuleType.Pathway)) {

				expregene = new ArrayList<>();
				expregene.addAll(orthologs);
				expregene.add(ec_number);
				expregene.add(reaction);
				ret = new HashSet<>();
			}
		}

		if(!expregene.isEmpty())
			ret = new HashSet<>();

		return ret;
	}


	/**
	 * @param orthologs
	 * @param reaction
	 * @param module
	 * @param moduleType
	 * @param mic
	 * @return
	 * @throws Exception
	 */
	private List<GeneAssociation> getdefinition(Set<String> orthologs, String reaction, String module, ModuleType moduleType, ModuleCI mic) throws Exception {

		//		List<String> orthologs = new ArrayList<>();
		//		{
		//			String query = "";
		//
		//			if(moduleType.equals(ModuleType.Complex)) {
		//
		//				boolean keep = false;
		//				for(String m : orthologsOfInterest) {
		//
		//					Pattern pattern = Pattern.compile("^ORTHOLOGY");
		//					Matcher matcher = pattern.matcher(m);
		//					if(matcher.find())
		//						keep = true;
		//
		//					pattern = Pattern.compile("^CLASS");
		//					matcher = pattern.matcher(m);
		//					if(matcher.find())
		//						keep = false;
		//
		//					if(keep)
		//						query = query.concat(m);
		//				}
		//			}
		//			else {
		//
		//				query = orthologsOfInterest[reaction2];
		//			}
		//
		//			Pattern pattern = Pattern.compile("K\\d{5}");
		//			Matcher matcher = pattern.matcher(query);
		//
		//			pattern = Pattern.compile("K\\d{5}");
		//			matcher = pattern.matcher(query);
		//			while(matcher.find()) {
		//
		//				orthologs.add(matcher.group());
		//			}
		//		}

		if(moduleType.equals(ModuleType.Complex)) {

			if(orthologs.size()==this.orthologsReactions.get(reaction).size()) { 

				orthologs.removeAll(this.orthologsReactions.get(reaction));

				if(orthologs.size()>0)
					return null;
			}
			else
				return null;
		} 
		else
			orthologs.retainAll(this.orthologsReactions.get(reaction));

		InputStream is = new ByteArrayInputStream(mic.getDefinition().getBytes());
		KEGGOrthologyParser parser = new KEGGOrthologyParser(is);

		List<List<String>> ret = parser.parseDefinition();
		
		logger.info("{}\t{}\t{} ", mic.getModule(),moduleType,ret);

		List<GeneAssociation> geneAssociationList = null;

		if(moduleType.equals(ModuleType.Pathway))
			geneAssociationList = this.getFinalRule(ret, reaction, mic, orthologs);
		else
			geneAssociationList = this.getFinalRule(0, ret, reaction, mic);

		return geneAssociationList;
	}

	/**
	 * @param moduleEntries
	 * @param module
	 * @param moduleType
	 * @return
	 */
	private ModuleCI parseModuleCI(KeggModulesParser k, String module, ModuleType moduleType) {

		ModuleCI mic = new ModuleCI(module, moduleType);

		mic.setDefinition(k.getDefinition());
		mic.setName(k.getName());

		return mic;
	}

	// Funcoes auxiliares
	/**
	 * @param frase
	 * @return
	 */
	private List<String> funcAux(String frase) {

		ArrayList<String> dividemais = new ArrayList<String>();
		dividemais.addAll(Arrays.asList(frase.split("\\+")));
		return dividemais;
	}

	/**
	 * @param listalista
	 * @return
	 */
	private List<List<String>> normalize(List<List<String>> lists){

		ArrayList<String> res = new ArrayList<String>();
		ArrayList<List<String>> result = new ArrayList<List<String>>();

		for(List<String> out : lists) {

			for (String s : out) {

				res.addAll(funcAux(s));
			}
			result.add(res);
			res = new ArrayList<String>();
		}
		return result;
	}

	/**
	 * @param definition
	 * @param reaction
	 * @param module
	 * @param moduleType
	 * @param orthologs
	 * @return
	 * @throws Exception
	 */
	private List<GeneAssociation> getFinalRule(List<List<String>> definition, String reaction, ModuleCI mic, Set<String> orthologs) throws Exception {

		List<List<String>> normalizedDefinition = this.normalize(definition);
		int index = -1;

		for(int i = 0; i<normalizedDefinition.size(); i++) {

			List<String> orthologsList = normalizedDefinition.get(i);

			//if(this.initialOrthologsReactions.get(reaction).containsAll(orthologsList))
			if(orthologsList.containsAll(orthologs))
				index = normalizedDefinition.indexOf(orthologsList);
		}

		return this.getFinalRule(index, definition, reaction, mic);
	}

	/**
	 * @param index
	 * @param definition
	 * @param reaction
	 * @param mic
	 * @return
	 * @throws Exception
	 */
	private List<GeneAssociation> getFinalRule(int index, List<List<String>> definition, String reaction, ModuleCI mic) throws Exception {

		List<GeneAssociation> gene_rules = new ArrayList<>();

		if(index>-1) {

			List<String> pathways = KeggAPI.getPathwaysIDByReaction(reaction);
			List<String> pathways_module = KeggAPI.getPathwaysByModule(mic.getModule());
			pathways.retainAll(pathways_module);
			mic.setPathways(pathways_module);

			String express = definition.get(index).toString().replaceAll(",", " or").replaceAll("\\+", " and ").replaceAll("\\[", "").replaceAll("\\]", "").trim();

			String[] rules = express.split("or");

			for(String rule : rules) {

				rule = rule.trim();

				GeneAssociation gene_rule = new GeneAssociation(mic);

				String[] and_rules = rule.split(" and ");
				
				for(String and_rule : and_rules)
					if(this.orthologsReactions.get(reaction).contains(and_rule.trim()))
						gene_rule.addGene(and_rule);

				if(gene_rule.getGenes().size()>0)
					gene_rules.add(gene_rule);
			}

			logger.info("GPR RULE:\t {}",gene_rules);
			return gene_rules;
		}
		return null;
	}

	/**
	 * @param orthologs
	 * @return
	 * @throws Exception 
	 */
	private List<String> getOrthologsWihtoutECnumber(List<String> orthologs) throws Exception {

		List<String> ret = new ArrayList<String>();

		for(String ortholog : orthologs) {

			//if(KeggAPI.getLinkECnumbersFromOrthology(ortholog).isEmpty())
			Set<String> s = KeggAPI.getECnumbersByOrthology(ortholog);

			if(s.isEmpty())
				ret.add(ortholog);
		}
		return ret;
	}
}