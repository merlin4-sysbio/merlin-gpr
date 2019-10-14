
package pt.uminho.ceb.biosystems.merlin.gpr.rules.core;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.ceb.biosystems.merlin.core.containers.gpr.ReactionsGPR_CI;
import pt.uminho.ceb.biosystems.merlin.dataAccess.InitDataAccess;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelReactionsServices;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.pair.Pair;

/**
 * @author ODias
 *
 */
public class FilterModelReactions {

	private static final Logger logger = LoggerFactory.getLogger(FilterModelReactions.class);

	private Map<String, List<String>> databaseEnzymesReactions;
	private Map<String, String> annotations; 
	private boolean isCompartimentalised;
	private Set<String> removed , kept, no_gpr;
	private Set<String> keptWithDifferentAnnotation;

	private String workspaceName;


	/**
	 * @param workspaceName
	 * @param originalReactions
	 */
	public FilterModelReactions(String workspaceName, boolean isCompartimentalised) { 
		this.workspaceName = workspaceName;
		this.databaseEnzymesReactions = new HashMap<>();
		this.isCompartimentalised = isCompartimentalised;

		this.removed = (new HashSet<>());
		this.kept = (new HashSet<>());
		this.keptWithDifferentAnnotation = (new HashSet<>());
		this.annotations  = new HashMap<>();
		this.no_gpr = (new HashSet<>());

	}


	/**
	 * @param keepReactionsWithNotes
	 * @param keepManualReactions
	 * @throws Exception
	 */
	public void removeReactionsFromModel(boolean keepReactionsWithNotes, boolean keepManualReactions) throws Exception {

		ModelReactionsServices.removeReactionsFromModelByBooleanRule(this.workspaceName, 
				new ArrayList<String>(this.removed), keepReactionsWithNotes, keepManualReactions);
	}


	/**
	 * @throws Exception
	 */
	public void setModelGPRsFromTool() throws Exception {
		
		ModelReactionsServices.updateBooleanRuleAndNotes(this.workspaceName, new ArrayList<String>(this.kept),
				this.annotations, "automatic GPR");
		
		ModelReactionsServices.updateBooleanRuleAndNotes(this.workspaceName, new ArrayList<String>(this.keptWithDifferentAnnotation), 
				this.annotations, "New Annotation. GPR set from tool");

	}
	

	public List<String> getReactionsFromModel(String databaseName, boolean isCompartimentalised) throws Exception {

		List<String> ret = new ArrayList<>();

		List<Pair<Integer, String>> result = InitDataAccess.getInstance().getDatabaseService(databaseName).getReactionHasEnzyme(isCompartimentalised);
		Pair<Integer, String> list;
		
		for(int i=0; i<result.size(); i++) {	
			
			list = result.get(i);
			
			List<String> reactions = new ArrayList<>();

			if(this.databaseEnzymesReactions.containsKey(String.valueOf(list.getB())));
				reactions = this.databaseEnzymesReactions.get(String.valueOf(list.getB()));

			reactions.add(String.valueOf(list.getA()));
			ret.add(String.valueOf(list.getA()));

			this.databaseEnzymesReactions.put(list.getB(), reactions);
		}
	
		return ret;
	}
	
	

	/**
	 * @param runGPRsAssignment
	 * @throws Exception 
	 */
	public void filterReactions(Map<String, ReactionsGPR_CI> runGPRsAssignment, boolean isCompartentalised ) throws Exception {
		
		Map<String, Set<String>> gpr_map = new HashMap<>();
		Set<String> ecs = new HashSet<>();

		for(String reaction : runGPRsAssignment.keySet()) {

			for(String key : runGPRsAssignment.get(reaction).getProteins().keySet()) {

				ecs.add(key);

				if(runGPRsAssignment.get(reaction).getProteins().get(key).isECnumberValid()) {

					this.annotations.put(reaction, runGPRsAssignment.get(reaction).getProteins().get(key).getGeneRule());

					Set<String> reactions = new HashSet<>();

					if(gpr_map.containsKey(key))
						reactions = gpr_map.get(key);

					reactions.add(reaction);

					gpr_map.put(key, reactions);
				}
			}
		}

		List<String> reactions = this.getReactionsFromModel(this.workspaceName,isCompartentalised);

		for(String ec : gpr_map.keySet()) {

			if(this.databaseEnzymesReactions.containsKey(ec)) {

				List<String> r = this.databaseEnzymesReactions.get(ec);
				List<String> r2 = new ArrayList<>();

				for(String react : r) {

					for(String gpr_react : gpr_map.get(ec)) {

						if(react.contains(gpr_react)) {

							r2.add(react);
							kept.add(react);
							this.annotations.put(react, this.annotations.get(gpr_react));
						}
						else {

							removed.add(react);
						}
					}
				}
			}
			else {

				for(String gpr_react : gpr_map.get(ec)) {

					for(String react : reactions) {

						if(react.contains(gpr_react)) {

							keptWithDifferentAnnotation.add(react);
							this.annotations.put(react, this.annotations.get(gpr_react));
						}
					}
				}
			}
		}

		for(String ec : this.databaseEnzymesReactions.keySet()) {

			if(!gpr_map.containsKey(ec)) {

				if(this.databaseEnzymesReactions.containsKey(ec)) {

					List<String> r = this.databaseEnzymesReactions.get(ec);

					for(String react : r) {

						this.no_gpr.add(react);
					}
				}
			}
		}

		this.keptWithDifferentAnnotation.removeAll(this.kept);
		this.removed.removeAll(this.kept);
		this.removed.removeAll(this.keptWithDifferentAnnotation);
		this.no_gpr.removeAll(this.kept);
		this.no_gpr.removeAll(this.keptWithDifferentAnnotation);
		this.removed.removeAll(this.no_gpr);

		String s = "Removed\t"+removed.size()+"\t"+removed+
				"\nKept\t"+kept.size()+"\t"+kept+
				"\nKept new annotation\t"+keptWithDifferentAnnotation.size()+"\t"+keptWithDifferentAnnotation+
				"\nNo GPR\t"+no_gpr.size()+"\t"+no_gpr;
		
		s+="\n\n";
		for(String r : this.annotations.keySet())
			s+="reaction:\t"+r+"\t"+this.annotations.get(r)+"\n";

		logger.info("Integration report: {}\nEnd Report.", s);
		
		
		
	}

}

