package pt.uminho.sysbio.merlin.gpr.rules.core.input;

public class KeggModulesParser extends AbstractKEGGFlatFileParser{

	public KeggModulesParser(String flatfile) {
		
		super(flatfile);
		super.parseContent();
	}
	
	
	public String getEntry() {
		
		String ret = null;
		
		int index = this.getTabIndex("entry");

		if(index>=0)
			ret = this.getContent(index);
		
		return ret;
	}
	
	public String getName() {
		
		String ret = null;
		
		int index = this.getTabIndex("name");
		
		if(index>=0)
			ret = this.getContent(index);
	
		return ret;
	}
	
	
	public String getDefinition() {
		
		String ret = null;
		
		int index = this.getTabIndex("definition");
		
		if(index>=0)
			ret = this.getContent(index);
	
		return ret;
	}

	public String getOrthology() {
		
		String ret = null;
		
		int index = this.getTabIndex("orthology");
		
		if(index>=0)
			ret = this.getContent(index);
	
		return ret;
	}

}
