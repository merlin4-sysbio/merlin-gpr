package pt.uminho.ceb.biosystems.merlin.gpr.rules.grammar;

import static org.junit.Assert.*;

import java.io.ByteArrayInputStream;

import org.junit.BeforeClass;
import org.junit.Test;

public class TestKEGGOrthologyParser {

	private static String[] examples = {"((((((((A))))))))",
										"A B",
										"(A, B, C) D",
										"(A, B) (C, D)",
										"A+B",
										"(A)+B",
										"(A,B)+C",
										"(A,B)+(C,D)",
										"(A,B)+(C,D) (X,Y)+(Z,W)",
										"(A,B)+((Z,W),(X,Y))",
										"(A,B)+((Z,W)+C,(X,Y)+D)",
										"(A,B)+((Z,W)+C,(X,Y,(X))+D)+M",
										"((A B), C)",
										"((A B), (C D))",
										"((A B)+Y, (C D)+Z)",
										"((A B)+Y, (C (O M G))+Z)+W",
										"((A B)+Y, (C (O M)-G)+Z)+W",
										"((A B C D,E,F,G),H)",
										"((A B,C D E,F O,M,G),H)",
										"((A B,(C D) E),H)",};
	
	private static String[] KEGGExamples = {"(K02777,K02817)+K02818+K02819", // M00270
											"K01655 K01681 K01705 K05824 (K00825,K00838) (K00143,K14085) K00293 K00290", // M00030
											"(((K01850,K04092,K14187,K04093,K04516,K06208,K06209,K13853)+(K01713,K04518,K05359)),K14170) (K00832,K00838)", // M00024
											"(((K01850,K04092,K14170,K04093,K04516,K06208,K06209,K13853)+(K00210,K04517)),K14187) (K00815,K00832,K00838)", // M00025
											"(K01655,K02594,K10977) ((K01681 K01705),K16792+K16793) (K05824,K10978)", // M00433
											"(K00169+K00170+K00171+K00172,K03737) ((K01007,K01006)+K01595,K01959+K01960,K01958) K00024 (K01676,K01679,K01677+K01678) (K00239+K00240-K00241-K00242,K00244+K00245+K00246+K00247) (K01902+K01903) (K00174+K00175-K00177-K00176) K00031 (K01681,K01682) (K15230+K15231,K15232+K15233+K15234)", //M00173
											"(K00169+K00170+K00171+K00172,K03737) K01007 K01595 K00024 (K01676,K01679,K01677+K01678) (K00239+K00240-K00241-K00242,K00244+K00245+K00246+K00247) (K01902+K01903) K15038 K14465 (K14466,K14467) K14534 (K15016,K01715+(K01782,K01825,K07516)) K00626" //M00374
											
										   };
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		KEGGOrthologyParser.VERBOSE = false;
	}

	@Test
	public void testExamples() {
		for (int i = 0; i < examples.length; i++) {
			ByteArrayInputStream stream = new ByteArrayInputStream(examples[i].getBytes());
			KEGGOrthologyParser koParser = new KEGGOrthologyParser(stream);
			try {
				System.out.println("PARSE RESULT: " + koParser.parseDefinition());
				
				System.out.println(examples[i] + " Passed !");
			} catch (ParseException pEx) {
				System.err.println("PEX: " + pEx.getMessage());
				System.out.println(examples[i] + " Failed !");
			}
			
		}
		
		assertEquals(true, true);
	}
	
	@Test
	public void testKEGGExamples() {
		for (int i = 0; i < KEGGExamples.length; i++) {
			ByteArrayInputStream stream = new ByteArrayInputStream(KEGGExamples[i].getBytes());
			KEGGOrthologyParser koParser = new KEGGOrthologyParser(stream);
			try {
				System.out.println("PARSE RESULT: " + koParser.parseDefinition());
				
				System.out.println(KEGGExamples[i] + " Passed !");
			} catch (ParseException pEx) {
				System.err.println("PEX: " + pEx.getMessage());
				System.out.println(KEGGExamples[i] + " Failed !");
			}
			
		}
		
		assertEquals(true, true);
	}

}
