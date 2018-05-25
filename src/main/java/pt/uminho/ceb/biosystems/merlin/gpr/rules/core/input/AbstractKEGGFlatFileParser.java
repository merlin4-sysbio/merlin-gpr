package pt.uminho.ceb.biosystems.merlin.gpr.rules.core.input;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public abstract class AbstractKEGGFlatFileParser {
	protected String flatfile_ = null;
	protected HashMap<Integer, String> tabName_ = null;
	protected HashMap<Integer, String> tabContent_ = null;
	private final String TAB_REGEX = "\n(\\w+)\\s+";
	private final String END_OF_FILE_DELIMITER = "///";
	
	public AbstractKEGGFlatFileParser( String flatfile) {
		this.flatfile_ = flatfile;
	}
	public String getFlatFile() {
		return this.flatfile_;
	}
	protected int getTabIndex( String name) {
		String dcName = name.toLowerCase();
		for ( Integer i : this.tabContent_.keySet()) {
			if ( this.tabName_.get(i).equals(dcName)) return i;
		}
		return -1;
	}
	protected void parseContent() {
		if ( flatfile_ == null) return;
		this.tabContent_ = new HashMap<Integer, String> ();
		this.tabName_ = new HashMap<Integer, String> ();
		//PARSE FIRST LINE
		String[] header = this.flatfile_.split("\\n")[0].split("\\s+");
		tabName_.put(0, header[0].toLowerCase());
		String concat = header[1];
		for(int i = 2; i<header.length; i++ ) {
			
			concat = concat.concat("\t"+header[i]);			
		}
		tabContent_.put(0, concat);
		//PARSE SECOND TAB TO LAST - 1 (MISSING LAST AND FIRST)

		Pattern tabPattern = Pattern.compile( TAB_REGEX);
		Matcher parser = tabPattern.matcher( this.flatfile_);
		int index = 1;
		int start = 0;
		int end = 0;
		String tabName = null;
		while ( parser.find() ) {
			end = parser.start();
			if ( start != 0) {
				tabName_.put( index, tabName.toLowerCase());
				tabContent_.put( index, this.flatfile_.substring( start, end));
				index++;
			}
			tabName = parser.group( 1);
			start = parser.end();
		}
		//PARSE LAST
		tabName_.put( index, tabName.toLowerCase());
		int endOfFile = this.flatfile_.indexOf( END_OF_FILE_DELIMITER, 0);
		tabContent_.put( index, this.flatfile_.substring( start, --endOfFile));
		//	System.out.println(tabContent_);
		//	System.out.println(tabName_);
	}
	public int tabCount() {
		return this.tabContent_.size();
	}
	public String getContent( int index) {
		return this.tabContent_.get( index);
	}
	public void readFile( String filepath) {
		try {
			FileInputStream fstream = new FileInputStream( filepath);
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String readLine;
			StringBuilder sb = new StringBuilder();
			while ((readLine = br.readLine()) != null) {
				sb.append( readLine).append("\n");
			}
			br.close();
			this.flatfile_ = sb.toString();
			this.parseContent();
		} catch ( FileNotFoundException fnfEx) {
		} catch ( IOException ioEx) {
		}
	}
}

