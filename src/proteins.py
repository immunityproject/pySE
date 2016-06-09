pdbs = {
	'RT': 'pdbs/RT_1DLO_N.pdb',
	'P24': 'pdbs/p24_Hexa_3H4E_N.pdb',
	'P17': 'pdbs/p17_2GOL_N.pdb',
	'GP120': 'pdbs/gp120_3JWD_N.pdb',
	'INT': 'pdbs/Integrase_1BIS_N.pdb',
	'NEF': 'pdbs/Nef_1EFN_N.pdb',
	'PRO': 'pdbs/protease_Dimer_3IXO_N.pdb',
	'REV': 'pdbs/Rev_3LPH_N.pdb',
	'TAT': 'pdbs/tat_3MI9_N.pdb',
}

all_aminos = {
	"P24":{'A':"PIVQNLQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRLHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTHNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQExxNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQG",'B':"PIVQNLQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRLHPxxxxxxxxGQMREPRGSDIAGTTSTLQEQIGWMTHNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAExxxxxxxNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQ",'C':"PIVQNLQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRLHPVxxxxxxxGQMREPRGSDIAGTTSTLQEQIGWMTHNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQ",'D':"PIVQNLQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRLHPxxxxxxxxGQMREPRGSDIAGTTSTLQEQIGWMTHNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQxxxxxxxxxxxTLLVQNANPDCKTILKALGPGATLEEMMTACQ",'E':"PIVQNLQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRLHxxxxxxxxxxxxREPRGSDIAGTTSTLQEQIGWMTHNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQ",'F':"PIVQNLQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRLHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTHNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQ"},
	"GP120":{'A':"TEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTGxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxGSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSSQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTGAxxxxxxxxxxxxxxxxxxxxxxxxxxxxxGHCNIARAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHWFNCGGEFFYCNSTQLFNSTWFNSTxxxxxxxxxEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKA",'B':"TEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTGxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxGSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSSQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTGAxxxxxxxxxxxxxxxxxxxxxxxxxxxxxGHCNIARAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHWFNCGGEFFYCNSTQLFNSTWFxxxxxxxxxxxxxGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPT",'C':"MKKVVLGKKGDTVELTCTASQKKSIQFHWKNSNQIKILGNQGSFLTKGPSKLNDRADSRRSLWDQGNFPLIIKNLKIEDSDTYICEVEDQKEEVQLLVFGLTANSDTHLLQGQSLTLTLESPPGSSPSVQCRSPRGKNIQGGKTLSVSQLELQDSGTWTCTVLQNQKKVEFKIDIVVLAFQKAS",'D':"KKVVLGKKGDTVELTCTASQKKSIQFHWKNSNQIKILGNQGSFLTKGPSKLNDRADSRRSLWDQGNFPLIIKNLKIEDSDTYICEVEDQKEEVQLLVFGLTANSDTHLLQGQSLTLTLESPPGSSPSVQCRSPRGKNIQGGKTLSVSQLELQDSGTWTCTVLQNQKKVEFKIDIVVLAFQKAS"},

	"INT":{'A':"CSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTVHTDNGSNFTSTTVKAACWWAGIKQEFGIPxxxxxxxxIESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQ",'B':"CSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTVHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGVIESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQ"},

	"NEF":{'A':"ALFVALYDYEAITEDDLSFHKGEKFQILNSSEGDWWEARSLTTGETGYIPSNYVAPV",'C':"ALFVALYDYEAITEDDLSFHKGEKFQILNSSEGDWWEARSLTTGETGYIPSNYVAPV",'B':"RPQVPLRPMTYKAAVDLSHFLKEKGGLEGLIHSQRRQDILDLWIYHTQGYFPDWQNYTPGPGVRYPLTFGWCYKLVPVxxxxxxxxxxxxxxxxxxxxxxxxxxxxxREVLEWRFDSRLAFHHVARELHPEYF",'D':"RPQVPLRPMTYKAAVDLSHFLKEKGGLEGLIHSQRRQDILDLWIYHTQGYFPDWQNYTPGPGVRYPLTFGWCYKLVPVxxxxxxxxxxxxxxxxxxxxxxxxxxxxxREVLEWRFDSRLAFHHVARELHPEYF"},

	"P17":{'A':"SVLSGGELDKWEKIRLRPGGKKQYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTIAVLYCVHQRIDVKDTKEALDKIEEE"},

	"PRO":{'A':"PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVGQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF",'B':"PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVGQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"},

	"REV":{'A':"DEDSLKAVRLIKFLYQSNPPPNPEGTRQARRNRRRRWRERQRQIHSISERIRSTYLGRSAEP",'B':"DEDSLKAVRLIKFLYQSNPPPNPEGTRQARRNRRRRWRERQRQIHSISERIRSTY",'C':"SDEDSLKAVRLIKFLYQSNPPPNPEGTRQARRNRRRRWRERQRQIHSISERIRSTYLG",'D':"SDEDSLKAVRLIKFLYQSNPPPNPEGTRQARRNRRRRWRERQRQIHSISERIRSTYL"},

	"RT":{'A':"PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKQKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQCSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYAGIKVRQLCKLLRGTKALTEVVPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMKGAHTNDVKQLTEAVQKIATESIVIWGKTPKFKLPIQKETWEAWWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIIGAETFYVDGAANRETKLGKAGYVTDRGRQKVVPLTDTTNQKTELQAIHLALQDSGLEVNIVTDSQYALGIIQAQPDKSESELVSQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDGLVSAGI",'B':"PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFKKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDxxxxxxxxxxxxGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYPGIKVRQLSKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWETWWTEYWQATWIPEWEFVNTPPLVKLWY"},
	"TAT":{'A':"VECPFCDEVSKYEKLAKIGQGTFGEVFKARHRKTGQKVALKKVLMENEKEGFPITALREIKILQLLKHENVVNLIEICRTKxxxxxxxxGSIYLVFDFCEHDLAGLLSNVLVKFTLSEIKRVMQMLLNGLYYIHRNKILHRDMKAANVLITRDGVLKLADFGLARAFSLAKNSQPNRYpNRVVTLWYRPPELLLGERDYGPPIDLWGAGCIMAEMWTRSPIMQGNTEQHQLALISQLCGSITPEVWPNVDNYELYEKLELVKGQKRKVKDRLKAYVRDPYALDLIDKLLVLDPAQRIDSDDALNHDFFWSDPMPSDLKGMLSTHLTSMFEYLAPPRR",'B':"NNNKRWYFTREQLENSPSRRFGVDPDKELSYRQQAANLLQDMGQRLNVSQLTINTAIVYMHRFYMIQSFTQFPGNSVAPAALFLAAKVEEQPKKLEHVIKVAHTCLHPQESLPDTRSEAYLQQVQDLVILESIILQTLGFELTIDHPHTHVVKCTQLVRASKDLAQTSYFMATNSLHLTTFSLQYTPPVVACVCIHLACKWSNWEIPVSTDGKHWWEYVDATVTLELLDELTHEFLQILEKTPNRLxxxxxxxxC",'C':"MEPVDPRLEPWKHPGSQPKTACTNCYCKKCCFHCQVCFITKALGISYGR"}
}

chains = {
	'P24': [['A','B','C','D','E','F']],
	'GP120': [['A','B'],['C','D']],
	'P17': [['A']],
	'RT': [['A','B']],
	'INT':[['A','B']],
	'NEF':[['A','C'],['B','D']],
	'PRO':[['A','B']],
	'REV':[['A','B','C','D']],
	'TAT':[['A'],['B'],['C']],
	'ZIKA_E': [['A','C','E']],
	'ZIKA_M': [['B','D','F']],
}

chain_bases = {

	'P24': {'A': 1,
	        'B': 1,
	        'C': 1,
	        'D': 1,
	        'E': 1,
	        'F': 1},
	'GP120': {'A': 31,
	          'B': 31,
	          'C': 1000,
	          'D': 1001},
	'P17': {'A': 6
	        },
	'RT': {'A': 1,
	       'B': 1},
	'INT': {'A': 56,
	        'B': 56},
	'NEF': {'A': 85,
	        'C': 85,
	        'B': 71,
	        'D': 71},
	'PRO': {'A': 1,
	        'B': 1},
	'REV': {'A': 9,
	        'B': 9,
	        'C': 9,
	        'D': 9},
	'TAT': {'A': 8,
	        'B': 7,
	        'C': 1},

}

codes = {"ALA": "A",
         "ARG": "R",
         "ASN": "N",
         "ASP": "D",
         "ASX": "B",
         "CYS": "C",
         "GLU": "E",
         "GLN": "Q",
         "GLX": "Z",
         "GLY": "G",
         "HIS": "H",
         "ILE": "I",
         "LEU": "L",
         "LYS": "K",
         "MET": "M",
         "PHE": "F",
         "PRO": "P",
         "SER": "S",
         "THR": "T",
         "TRP": "W",
         "TYR": "Y",
         "VAL": "V",
         "TPO": "X"}

#manual insert hack
# 'ZIKA_M': {'VI10': {'chains': ['A', 'B'], 'reference_value': 2.33, 'last': 505, 'conf': 1.0, 'first': 496},

proteins = {'ZIKA_E':{},'RT': {'VI10': {'chains': ['A', 'B'], 'reference_value': 2.33, 'last': 505, 'conf': 1.0, 'first': 496}, 'ER11': {'chains': ['A', 'B'], 'reference_value': 2.17, 'last': 448, 'conf': 1.0, 'first': 438}, 'IY10': {'chains': ['A', 'B'], 'reference_value': 2.2, 'last': 318, 'conf': 0.8, 'first': 309}, 'GA10': {'chains': ['A', 'B'], 'reference_value': 2.04, 'last': 445, 'conf': 1.0, 'first': 436}, 'VY10': {'chains': ['A', 'B'], 'reference_value': 1.84, 'last': 127, 'conf': 1.0, 'first': 118}, 'TY9': {'chains': ['A', 'B'], 'reference_value': 1.72, 'last': 115, 'conf': 1.0, 'first': 107}, 'KK10': {'chains': ['A', 'B'], 'reference_value': 1.71, 'last': 82, 'conf': 1.0, 'first': 73}, 'IV9': {'chains': ['A', 'B'], 'reference_value': 2.19, 'last': 317, 'conf': 1.0, 'first': 309}, 'GK9': {'chains': ['A', 'B'], 'reference_value': 1.3, 'last': 101, 'conf': 1.0, 'first': 93}, 'QR9': {'chains': ['A', 'B'], 'reference_value': 2.35, 'last': 277, 'conf': 1.0, 'first': 269}, 'GI9': {'chains': ['A', 'B'], 'reference_value': 1.57, 'last': 341, 'conf': 1.0, 'first': 333}, 'IK10': {'chains': ['A', 'B'], 'reference_value': 2.07, 'last': 350, 'conf': 1.0, 'first': 341}, 'EY9': {'chains': ['A', 'B'], 'reference_value': 2.06, 'last': 457, 'conf': 1.0, 'first': 449}, 'DI8': {'chains': ['A', 'B'], 'reference_value': 2.38, 'last': 505, 'conf': 1.0, 'first': 498}, 'YI9': {'chains': ['A', 'B'], 'reference_value': 1.45, 'last': 135, 'conf': 1.0, 'first': 127}, 'AK11': {'chains': ['A', 'B'], 'reference_value': 2.28, 'last': 43, 'conf': 1.0, 'first': 33}, 'VL9': {'chains': ['A', 'B'], 'reference_value': 1.93, 'last': 325, 'conf': 1.0, 'first': 317}, 'PW10': {'chains': ['A', 'B'], 'reference_value': 1.93, 'last': 401, 'conf': 1.0, 'first': 392}, 'KY9': {'chains': ['A', 'B'], 'reference_value': 2.25, 'last': 181, 'conf': 1.0, 'first': 173}, 'TI8': {'chains': ['A', 'B'], 'reference_value': 1.43, 'last': 135, 'conf': 1.0, 'first': 128}, 'IW9': {'chains': ['A', 'B'], 'reference_value': 1.54, 'last': 383, 'conf': 0.8888888888888888, 'first': 375}, 'AM9': {'chains': ['A', 'B'], 'reference_value': 2.19, 'last': 41, 'conf': 1.0, 'first': 33}, 'AK9': {'chains': ['A', 'B'], 'reference_value': 2.02, 'last': 166, 'conf': 1.0, 'first': 158}, 'NY9': {'chains': ['A', 'B'], 'reference_value': 2.1, 'last': 183, 'conf': 0.8888888888888888, 'first': 175}, 'LY12': {'chains': ['A', 'B'], 'reference_value': 1.78, 'last': 271, 'conf': 1.0, 'first': 260}, 'SM9': {'chains': ['A', 'B'], 'reference_value': 1.56, 'last': 164, 'conf': 1.0, 'first': 156}, 'EI9': {'chains': ['A', 'B'], 'reference_value': 2.18, 'last': 50, 'conf': 1.0, 'first': 42}, 'QK9': {'chains': ['A', 'B'], 'reference_value': 2.43, 'last': 528, 'conf': 1.0, 'first': 520}, 'YL9': {'chains': ['A', 'B'], 'reference_value': 2.2, 'last': 279, 'conf': 1.0, 'first': 271}, 'GL9': {'chains': ['A', 'B'], 'reference_value': 1.67, 'last': 26, 'conf': 1.0, 'first': 18}, 'DV9': {'chains': ['A', 'B'], 'reference_value': 1.88, 'last': 372, 'conf': 0.8888888888888888, 'first': 364}, 'KY11': {'chains': ['A', 'B'], 'reference_value': 2.24, 'last': 183, 'conf': 1.0, 'first': 173}, 'HY9': {'chains': ['A', 'B'], 'reference_value': 2.1, 'last': 183, 'conf': 0.8888888888888888, 'first': 175}, 'IL8': {'chains': ['A', 'B'], 'reference_value': 2.24, 'last': 12, 'conf': 1.0, 'first': 5}, 'IL9': {'chains': ['A', 'B'], 'reference_value': 2.26, 'last': 503, 'conf': 1.0, 'first': 495}}, 'TAT': {'IR11': {'chains': ['C'], 'reference_value': 2.29, 'last': 49, 'conf': 0.9090909090909091, 'first': 39}, 'FY10': {'chains': ['C'], 'reference_value': 2.12, 'last': 47, 'conf': 0.8, 'first': 38}, 'CC8': {'chains': ['C'], 'reference_value': 1.58, 'last': 37, 'conf': 1.0, 'first': 30}, 'EW10': {'chains': ['C'], 'reference_value': 2.68, 'last': 11, 'conf': 1.0, 'first': 2}}, 'INT': {'FY10': {'chains': ['A', 'B'], 'reference_value': 2.1, 'last': 194, 'conf': 1.0, 'first': 185}, 'TL9': {'chains': ['A', 'B'], 'reference_value': 1.83, 'last': 74, 'conf': 0.8888888888888888, 'first': 66}, 'KF9': {'chains': ['A', 'B'], 'reference_value': 1.68, 'last': 181, 'conf': 1.0, 'first': 173}, 'KY9': {'chains': ['A', 'B'], 'reference_value': 2.06, 'last': 194, 'conf': 1.0, 'first': 186}, 'AK10': {'chains': ['A', 'B'], 'reference_value': 1.92, 'last': 188, 'conf': 1.0, 'first': 179}, 'IY9': {'chains': ['A', 'B'], 'reference_value': 2.35, 'last': 143, 'conf': 0.8888888888888888, 'first': 135}}, 'PRO': {'RI8': {'chains': ['A', 'B'], 'reference_value': 1.72, 'last': 15, 'conf': 1.0, 'first': 8}, 'LI9': {'chains': ['A', 'B'], 'reference_value': 1.64, 'last': 84, 'conf': 1.0, 'first': 76}, 'RI10': {'chains': ['A', 'B'], 'reference_value': 1.89, 'last': 66, 'conf': 0.9, 'first': 57}, 'GL9': {'chains': ['A', 'B'], 'reference_value': 1.75, 'last': 76, 'conf': 0.8888888888888888, 'first': 68}, 'EW9': {'chains': ['A', 'B'], 'reference_value': 2.13, 'last': 42, 'conf': 1.0, 'first': 34}, 'TL11': {'chains': ['A', 'B'], 'reference_value': 1.59, 'last': 90, 'conf': 0.9090909090909091, 'first': 80}, 'DL9': {'chains': ['A', 'B'], 'reference_value': 1.83, 'last': 38, 'conf': 1.0, 'first': 30}, 'IV9': {'chains': ['A', 'B'], 'reference_value': 2.01, 'last': 11, 'conf': 1.0, 'first': 3}, 'KV8': {'chains': ['A', 'B'], 'reference_value': 1.52, 'last': 77, 'conf': 1.0, 'first': 70}}, 'REV': {'IY9': {'chains': ['A', 'B', 'C', 'D'], 'reference_value': 1.9, 'last': 63, 'conf': 0.8888888888888888, 'first': 55}, 'ER10': {'chains': ['A', 'B', 'C', 'D'], 'reference_value': 2.05, 'last': 66, 'conf': 0.9, 'first': 57}, 'KY10': {'chains': ['A', 'B', 'C', 'D'], 'reference_value': 1.16, 'last': 23, 'conf': 1.0, 'first': 14}}, 'GP120': {'VL11': {'chains': ['A', 'B'], 'reference_value': 2.44, 'last': 52, 'conf': 1.0, 'first': 42}, 'LY10': {'chains': ['A', 'B'], 'reference_value': 2.26, 'last': 61, 'conf': 1.0, 'first': 52}, 'SF9': {'chains': ['A', 'B'], 'reference_value': 1.19, 'last': 383, 'conf': 0.8888888888888888, 'first': 375}, 'TK10': {'chains': ['A', 'B'], 'reference_value': 2.63, 'last': 46, 'conf': 1.0, 'first': 37}, 'AY10': {'chains': ['A', 'B'], 'reference_value': 2.93, 'last': 40, 'conf': 0.8, 'first': 31}, 'MW9': {'chains': ['A', 'B'], 'reference_value': 1.83, 'last': 112, 'conf': 1.0, 'first': 104}, 'SY10': {'chains': ['A', 'B'], 'reference_value': 1.11, 'last': 384, 'conf': 0.9, 'first': 375}, 'KW11': {'chains': ['A', 'B'], 'reference_value': 2.24, 'last': 69, 'conf': 0.9090909090909091, 'first': 59}, 'VT10': {'chains': ['A', 'B'], 'reference_value': 2.6, 'last': 51, 'conf': 1.0, 'first': 42}, 'YW9': {'chains': ['A', 'B'], 'reference_value': 2.1, 'last': 69, 'conf': 0.8888888888888888, 'first': 61}, 'LI9': {'chains': ['A', 'B'], 'reference_value': 1.96, 'last': 424, 'conf': 1.0, 'first': 416}, 'AY9': {'chains': ['A', 'B'], 'reference_value': 2.94, 'last': 39, 'conf': 0.7777777777777778, 'first': 31}, 'SK9': {'chains': ['A', 'B'], 'reference_value': 2.5, 'last': 207, 'conf': 1.0, 'first': 199}, 'RW9': {'chains': ['A', 'B'], 'reference_value': 2.14, 'last': 427, 'conf': 1.0, 'first': 419}, 'SY9': {'chains': ['A', 'B'], 'reference_value': 1.97, 'last': 217, 'conf': 1.0, 'first': 209}, 'FY9': {'chains': ['A', 'B'], 'reference_value': 1.14, 'last': 384, 'conf': 1.0, 'first': 376}, 'DL9': {'chains': ['A', 'B'], 'reference_value': 2.66, 'last': 86, 'conf': 1.0, 'first': 78}}, 'NEF': {'LV10': {'chains': ['B', 'D'], 'reference_value': 1.4, 'last': 146, 'conf': 0.9, 'first': 137}, 'HQ10': {'chains': ['B', 'D'], 'reference_value': 1.38, 'last': 125, 'conf': 1.0, 'first': 116}, 'AL9': {'chains': ['B', 'D'], 'reference_value': 1.7, 'last': 91, 'conf': 1.0, 'first': 83}, 'YY9': {'chains': ['B', 'D'], 'reference_value': 1.73, 'last': 135, 'conf': 0.8888888888888888, 'first': 127}, 'RL9': {'chains': ['B', 'D'], 'reference_value': 1.75, 'last': 85, 'conf': 0.8888888888888888, 'first': 77}, 'KL10': {'chains': ['B', 'D'], 'reference_value': 1.81, 'last': 91, 'conf': 1.0, 'first': 82}, 'TL10': {'chains': ['B', 'D'], 'reference_value': 1.65, 'last': 137, 'conf': 1.0, 'first': 128}, 'RW8': {'chains': ['B', 'D'], 'reference_value': 1.23, 'last': 141, 'conf': 1.0, 'first': 134}, 'YT9': {'chains': ['B', 'D'], 'reference_value': 1.33, 'last': 128, 'conf': 1.0, 'first': 120}, 'WH10': {'chains': ['B', 'D'], 'reference_value': 1.92, 'last': 192, 'conf': 1.0, 'first': 183}, 'KL9': {'chains': ['B', 'D'], 'reference_value': 2.25, 'last': 100, 'conf': 1.0, 'first': 92}, 'RY10': {'chains': ['B', 'D'], 'reference_value': 1.89, 'last': 115, 'conf': 1.0, 'first': 106}, 'RY11': {'chains': ['B', 'D'], 'reference_value': 1.96, 'last': 115, 'conf': 0.9090909090909091, 'first': 105}, 'LL9': {'chains': ['B', 'D'], 'reference_value': 1.29, 'last': 145, 'conf': 0.8888888888888888, 'first': 137}, 'RI9': {'chains': ['B', 'D'], 'reference_value': 1.93, 'last': 114, 'conf': 1.0, 'first': 106}, 'PK8': {'chains': ['B', 'D'], 'reference_value': 1.74, 'last': 82, 'conf': 1.0, 'first': 75}, 'RM9': {'chains': ['B', 'D'], 'reference_value': 1.66, 'last': 79, 'conf': 1.0, 'first': 71}, 'RI10': {'chains': ['B', 'D'], 'reference_value': 2.01, 'last': 114, 'conf': 1.0, 'first': 105}, 'HW9': {'chains': ['B', 'D'], 'reference_value': 1.44, 'last': 124, 'conf': 1.0, 'first': 116}, 'RV9': {'chains': ['B', 'D'], 'reference_value': 1.75, 'last': 85, 'conf': 1.0, 'first': 77}, 'YF9': {'chains': ['B', 'D'], 'reference_value': 1.21, 'last': 143, 'conf': 0.8888888888888888, 'first': 135}, 'FL8': {'chains': ['B', 'D'], 'reference_value': 2.22, 'last': 97, 'conf': 1.0, 'first': 90}, 'AK9': {'chains': ['B', 'D'], 'reference_value': 1.78, 'last': 92, 'conf': 1.0, 'first': 84}, 'PL10': {'chains': ['B', 'D'], 'reference_value': 1.17, 'last': 145, 'conf': 1.0, 'first': 136}, 'VL10': {'chains': ['B', 'D'], 'reference_value': 1.83, 'last': 189, 'conf': 1.0, 'first': 180}, 'TM9': {'chains': ['B', 'D'], 'reference_value': 1.66, 'last': 79, 'conf': 0.8888888888888888, 'first': 71}, 'VY8': {'chains': ['B', 'D'], 'reference_value': 1.61, 'last': 81, 'conf': 1.0, 'first': 74}, 'TY11': {'chains': ['B', 'D'], 'reference_value': 1.29, 'last': 127, 'conf': 1.0, 'first': 117}, 'GL9': {'chains': ['B', 'D'], 'reference_value': 1.7, 'last': 91, 'conf': 0.6666666666666666, 'first': 83}, 'WF9': {'chains': ['B', 'D'], 'reference_value': 1.83, 'last': 191, 'conf': 1.0, 'first': 183}, 'TW9': {'chains': ['B', 'D'], 'reference_value': 1.37, 'last': 141, 'conf': 0.8888888888888888, 'first': 133}, 'KY11': {'chains': ['B', 'D'], 'reference_value': 1.96, 'last': 115, 'conf': 0.7272727272727273, 'first': 105}, 'QK10': {'chains': ['B', 'D'], 'reference_value': 1.84, 'last': 82, 'conf': 1.0, 'first': 73}}, 'P24': {'PY9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.52, 'last': 130, 'conf': 0.8888888888888888, 'first': 122}, 'KI8': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.84, 'last': 37, 'conf': 1.0, 'first': 30}, 'AV9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.89, 'last': 86, 'conf': 0.8888888888888888, 'first': 78}, 'FF9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.38, 'last': 40, 'conf': 1.0, 'first': 32}, 'AW11': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.85, 'last': 184, 'conf': 0.9090909090909091, 'first': 174}, 'QW9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.87, 'last': 184, 'conf': 1.0, 'first': 176}, 'FK10': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.09, 'last': 170, 'conf': 1.0, 'first': 161}, 'VI9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.26, 'last': 150, 'conf': 0.8888888888888888, 'first': 142}, 'KA9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.71, 'last': 78, 'conf': 1.0, 'first': 70}, 'QW11': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.27, 'last': 23, 'conf': 1.0, 'first': 13}, 'VV9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.41, 'last': 11, 'conf': 1.0, 'first': 3}, 'GL9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.78, 'last': 69, 'conf': 1.0, 'first': 61}, 'RL11': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.16, 'last': 172, 'conf': 0.9090909090909091, 'first': 162}, 'TL9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.89, 'last': 56, 'conf': 0.6666666666666666, 'first': 48}, 'GI11': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.15, 'last': 104, 'conf': 1.0, 'first': 94}, 'HA9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 2.59, 'last': 92, 'conf': 1.0, 'first': 84}, 'EL9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.55, 'last': 43, 'conf': 1.0, 'first': 35}, 'IW9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.03, 'last': 23, 'conf': 1.0, 'first': 15}, 'SV9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.12, 'last': 24, 'conf': 1.0, 'first': 16}, 'DL9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.61, 'last': 205, 'conf': 1.0, 'first': 197}, 'EW10': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.59, 'last': 80, 'conf': 1.0, 'first': 71}, 'HL9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.48, 'last': 20, 'conf': 1.0, 'first': 12}, 'RI8': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.35, 'last': 150, 'conf': 1.0, 'first': 143}, 'VL8': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.62, 'last': 43, 'conf': 1.0, 'first': 36}, 'KK10': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.16, 'last': 140, 'conf': 1.0, 'first': 131}, 'RV9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.66, 'last': 181, 'conf': 1.0, 'first': 173}, 'DA9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.03, 'last': 174, 'conf': 1.0, 'first': 166}, 'GY9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.82, 'last': 145, 'conf': 1.0, 'first': 137}, 'SL9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.96, 'last': 52, 'conf': 1.0, 'first': 44}, 'TW10': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.08, 'last': 117, 'conf': 1.0, 'first': 108}, 'RK10': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.16, 'last': 140, 'conf': 0.7, 'first': 131}, 'AM12': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.84, 'last': 185, 'conf': 1.0, 'first': 174}, 'EI8': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.87, 'last': 135, 'conf': 1.0, 'first': 128}, 'YL9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.97, 'last': 172, 'conf': 0.8888888888888888, 'first': 164}, 'EV9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.18, 'last': 36, 'conf': 1.0, 'first': 28}, 'VF9': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 1.65, 'last': 32, 'conf': 0.8888888888888888, 'first': 24}, 'KF11': {'chains': ['A', 'B', 'C', 'D', 'E', 'F'], 'reference_value': 0.71, 'last': 40, 'conf': 1.0, 'first': 30}}, 'P17': {'RY10': {'chains': ['A'], 'reference_value': 2.73, 'last': 29, 'conf': 0.9, 'first': 20}, 'GK9': {'chains': ['A'], 'reference_value': 2.77, 'last': 32, 'conf': 0.8888888888888888, 'first': 24}, 'LL8': {'chains': ['A'], 'reference_value': 2.04, 'last': 85, 'conf': 0.75, 'first': 78}, 'KK9': {'chains': ['A'], 'reference_value': 2.74, 'last': 26, 'conf': 1.0, 'first': 18}, 'TI9': {'chains': ['A'], 'reference_value': 2.24, 'last': 92, 'conf': 0.8888888888888888, 'first': 84}, 'HL9': {'chains': ['A'], 'reference_value': 2.5, 'last': 41, 'conf': 0.8888888888888888, 'first': 33}, 'RY11': {'chains': ['A'], 'reference_value': 2.18, 'last': 86, 'conf': 0.8181818181818182, 'first': 76}, 'IL10': {'chains': ['A'], 'reference_value': 2.48, 'last': 101, 'conf': 0.8, 'first': 92}, 'RK9': {'chains': ['A'], 'reference_value': 2.84, 'last': 28, 'conf': 0.8888888888888888, 'first': 20}, 'LY9': {'chains': ['A'], 'reference_value': 2.03, 'last': 86, 'conf': 0.7777777777777778, 'first': 78}, 'GI9': {'chains': ['A'], 'reference_value': 2.56, 'last': 19, 'conf': 0.8888888888888888, 'first': 11}, 'LF11': {'chains': ['A'], 'reference_value': 2.56, 'last': 44, 'conf': 0.9090909090909091, 'first': 34}, 'GL8': {'chains': ['A'], 'reference_value': 2.75, 'last': 31, 'conf': 0.875, 'first': 24}, 'WF9': {'chains': ['A'], 'reference_value': 2.56, 'last': 44, 'conf': 1.0, 'first': 36}, 'KW9': {'chains': ['A'], 'reference_value': 2.67, 'last': 36, 'conf': 0.8888888888888888, 'first': 28}, 'IK9': {'chains': ['A'], 'reference_value': 2.74, 'last': 27, 'conf': 1.0, 'first': 19}, 'TK8': {'chains': ['A'], 'reference_value': 2.19, 'last': 91, 'conf': 0.75, 'first': 84}, 'SL9': {'chains': ['A'], 'reference_value': 2.12, 'last': 85, 'conf': 0.7777777777777778, 'first': 77}}}




#proteins = {
#	'rt': {
#		'gk9-rt': {'first': 93, 'last': 101, 'reference_value': 2.77},
#		'iy10': {'first': 309, 'last': 318, 'reference_value': 2.20},
#		'ky11': {'first': 173, 'last': 183, 'reference_value': 2.24},
#		'ly12': {'first': 260, 'last': 271, 'reference_value': 1.78},
#		'iv9': {'first': 309, 'last': 317, 'reference_value': 2.01},
#		'ga10': {'first': 436, 'last': 445, 'reference_value': None},
#	},
#	'p24': {
#		'kf11': {'first': 30, 'last': 40, 'reference_value': 0.71},
#		'qw11': {'first': 13, 'last': 23, 'reference_value': 1.27},
#		'rv9': {'first': 173, 'last': 181, 'reference_value': None},
#		'ki8': {'first': 30, 'last': 37, 'reference_value': 0.84},
#		'ff9': {'first': 32, 'last': 40, 'reference_value': 0.38},
#		'vv9': {'first': 3, 'last': 11, 'reference_value': 1.41},
#		'hl9': {'first': 12, 'last': 20, 'reference_value': 1.48},
#		'iw9': {'first': 15, 'last': 23, 'reference_value': 1.03},
#		'sv9': {'first': 16, 'last': 24, 'reference_value': 1.12},
#		'vf9': {'first': 24, 'last': 32, 'reference_value': 1.65},
#		'ev9': {'first': 28, 'last': 36, 'reference_value': 1.18},
#		'el9': {'first': 35, 'last': 43, 'reference_value': 0.55},
#		'vl8': {'first': 36, 'last': 43, 'reference_value': 0.62},
#		'sl9': {'first': 44, 'last': 52, 'reference_value': 0.96},
#		'tl9': {'first': 48, 'last': 56, 'reference_value': 0.89},
#	},
#	'p17': {
#		'ry11': {'first': 76, 'last': 86, 'reference_value': 2.18},
#		'sl9': {'first': 77, 'last': 85, 'reference_value': 2.12},
#		'll8': {'first': 78, 'last': 85, 'reference_value': 2.04},
#		'ly9': {'first': 78, 'last': 86, 'reference_value': 2.03},
#		'tk8': {'first': 84, 'last': 91, 'reference_value': 2.19},
#		'ti9': {'first': 84, 'last': 92, 'reference_value': 2.24},
#		'il10': {'first': 92, 'last': 101, 'reference_value': 2.48},
#		'gi9': {'first': 11, 'last': 19, 'reference_value': 2.56},
#		'kk9': {'first': 18, 'last': 26, 'reference_value': 2.74},
#		'ik9': {'first': 19, 'last': 27, 'reference_value': 2.74},
#		'rk9': {'first': 20, 'last': 28, 'reference_value': 2.84},
#		'ry10': {'first': 20, 'last': 29, 'reference_value': 2.73},
#		'gl8': {'first': 24, 'last': 31, 'reference_value': 2.75},
#		'gk9': {'first': 24, 'last': 32, 'reference_value': 2.77},
#		'kw9': {'first': 28, 'last': 36, 'reference_value': 2.67},
#		'hl9': {'first': 33, 'last': 41, 'reference_value': 2.50},
#		'lf11': {'first': 34, 'last': 44, 'reference_value': 2.56},
#		'wf9': {'first': 36, 'last': 44, 'reference_value': 2.56},
#	}
#}
