/**
 * 
 */
package pt.uminho.sysbio.merlin.gpr.rules.core;

import java.util.HashSet;
import java.util.Set;

/**
 * @author ODias
 *
 */
public class ModuleConfirmation {

	private String module;
	private boolean hasECnumber;
	private boolean hasReaction;
	private Set<String> orthologs;

	/**
	 * 
	 */
	public ModuleConfirmation(String module) {
		
		this.setModule(module);
	}
	
	/**
	 * @param hasECnumber
	 * @param hasReaction
	 * @param orthologs
	 */
	public ModuleConfirmation(String module, boolean hasECnumber, boolean hasReaction, Set<String> orthologs) {
		
		this.setModule(module);
		this.setHasECnumber(hasECnumber);
		this.setHasReaction(hasReaction);
		this.setOrthologs(orthologs);
	}
	
	/**
	 * @param ortholog
	 */
	public void addOrtholog(String ortholog) {
		
		if(this.orthologs==null)
			this.orthologs = new HashSet<String>();
		
		this.orthologs.add(ortholog);
	}

	/**
	 * @return the hasECnumber
	 */
	public boolean isHasECnumber() {
		return hasECnumber;
	}

	/**
	 * @param hasECnumber the hasECnumber to set
	 */
	public void setHasECnumber(boolean hasECnumber) {
		this.hasECnumber = hasECnumber;
	}

	/**
	 * @return the hasReaction
	 */
	public boolean isHasReaction() {
		return hasReaction;
	}

	/**
	 * @param hasReaction the hasReaction to set
	 */
	public void setHasReaction(boolean hasReaction) {
		this.hasReaction = hasReaction;
	}

	/**
	 * @return the orthologs
	 */
	public Set<String> getOrthologs() {
		return orthologs;
	}

	/**
	 * @param orthologs the orthologs to set
	 */
	public void setOrthologs(Set<String> orthologs) {
		this.orthologs = orthologs;
	}

	/**
	 * @return the module
	 */
	public String getModule() {
		return module;
	}

	/**
	 * @param module the module to set
	 */
	public void setModule(String module) {
		this.module = module;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "ModuleConfirmation ["
				+ (module != null ? "module=" + module + ", " : "")
				+ "hasECnumber=" + hasECnumber + ", hasReaction=" + hasReaction
				+ ", " + (orthologs != null ? "orthologs=" + orthologs : "")
				+ "]";
	}

}
