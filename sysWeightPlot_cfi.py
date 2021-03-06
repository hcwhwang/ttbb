baseWeight = "tri*puweight*weight*mueffweight*eleffweight"
baseWeight2 = "tri*puweight*weight*mueffweight*eleffweight*csvweights2[0]"
baseWeightNOPU = "tri*weight*mueffweight*eleffweight"
mceventweightNOPU=[
#{"name":"CEN",       "tree":"nom2",    "var":"("+baseWeight+")"},
#{"name":"NoPU",      "tree":"nom2",    "var":"("+baseWeightNOPU+")"},
{"name":"csvweight", "tree":"nom2",    "var":baseWeightNOPU},
#{"name":"PU_Up",     "tree":"nom2",    "var":"(tri*weight*puweightUp*mueffweight*eleffweight)"},
#{"name":"PU_Down",   "tree":"nom2",    "var":"(tri*weight*puweightDown*mueffweight*eleffweight)"},
{"name":"JER_Up",    "tree":"nomJER_up",   "var":baseWeightNOPU},
#{"name":"JER_Down",  "tree":"nomJER_dw",   "var":baseWeightNOPU},
{"name":"JES_Up",    "tree":"nomJES_up",   "var":baseWeightNOPU},        
#{"name":"JES_Down",  "tree":"nomJES_dw",   "var":baseWeightNOPU},      
{"name":"MuEff_Up",     "tree":"nom2",    "var":"(tri*weight*mueffweight_up*eleffweight)"},
#{"name":"MuEff_Down",   "tree":"nom2",    "var":"(tri*weight*mueffweight_dn*eleffweight)"},
{"name":"ElEff_Up",     "tree":"nom2",    "var":"(tri*weight*mueffweight*eleffweight_up)"},
#{"name":"ElEff_Down",   "tree":"nom2",    "var":"(tri*weight*mueffweight*eleffweight_dn)"},
{"name":"Mu_Scale_Up",     "tree":"nomMu_up",    "var":baseWeightNOPU},
#{"name":"Mu_Scale_Down",     "tree":"nomMu_dw",  "var":baseWeightNOPU},
{"name":"El_Scale_Up",     "tree":"nomEl_up",    "var":baseWeightNOPU},
#{"name":"El_Scale_Down",     "tree":"nomEl_dw",  "var":baseWeightNOPU},
{"name":"Trigger_Up",   "tree":"nom2",        "var": "(weight*mueffweight*eleffweight*tri_up)"},
#{"name":"Trigger_Down", "tree":"nom2",        "var": "(weight*mueffweight*eleffweight*tri_dn)"},
{"name":"topPt_Up", "tree":"nom2",      "var":"("+baseWeightNOPU+"*topPtWeight)"},
#{"name":"topPt_Down", "tree":"nom2",    "var":"("+baseWeightNOPU+")"},
#{"name":"FSR_Down", "tree":"nom2",      "var":"("+baseWeightNOPU+")"},
{"name":"FSR_Up", "tree":"nom2",        "var":"("+baseWeightNOPU+")"},
#{"name":"ISR_Down", "tree":"nom2",      "var":"("+baseWeightNOPU+")"},
{"name":"ISR_Up", "tree":"nom2",        "var":"("+baseWeightNOPU+")"},
#{"name":"UE_Down", "tree":"nom2",       "var":"("+baseWeightNOPU+")"},
{"name":"UE_Up", "tree":"nom2",         "var":"("+baseWeightNOPU+")"},
#{"name":"hdamp_Down", "tree":"nom2",    "var":"("+baseWeightNOPU+")"},
{"name":"hdamp_Up", "tree":"nom2",      "var":"("+baseWeightNOPU+")"},
#{"name":"qcderdon_Down", "tree":"nom2", "var":"("+baseWeightNOPU+")"},
{"name":"qcderdon_Up", "tree":"nom2",   "var":"("+baseWeightNOPU+")"},
#{"name":"erdon_Down", "tree":"nom2",    "var":"("+baseWeightNOPU+")"},
{"name":"erdon_Up", "tree":"nom2",      "var":"("+baseWeightNOPU+")"},
#{"name":"gluonmove_Down", "tree":"nom2","var":"("+baseWeightNOPU+")"},
{"name":"gluonmove_Up", "tree":"nom2",  "var":"("+baseWeightNOPU+")"},
{"name":"MuF_Up", "tree":"nom2",        "var":"("+baseWeightNOPU+"*scaleWeightsUp[0])"},
#{"name":"MuF_Down", "tree":"nom2",      "var":"("+baseWeightNOPU+"*scaleWeightsDown[0])"},
#{"name":"MuR_Down", "tree":"nom2",      "var":"("+baseWeightNOPU+"*scaleWeightsDown[1])"},
{"name":"MuR_Up", "tree":"nom2",        "var":"("+baseWeightNOPU+"*scaleWeightsUp[1])"},
{"name":"PDF_Up", "tree":"nom2",        "var":"("+baseWeightNOPU+"*pdfWeights[0])"},
#{"name":"PDF_Down", "tree":"nom2",      "var":"("+baseWeightNOPU+"*pdfWeights[1])"},
{"name":"PDFAlphaS_Up", "tree":"nom2",  "var":"("+baseWeightNOPU+"*pdfWeights[2])"},
#{"name":"PDFAlphaS_Down", "tree":"nom2","var":"("+baseWeightNOPU+"*pdfWeights[3])"},
]
mceventweight=[
#{"name":"CEN",       "tree":"nom2",    "var":"("+baseWeight+")"},
#{"name":"NoPU",      "tree":"nom2",    "var":"("+baseWeightNOPU+")"},
{"name":"csvweight", "tree":"nom2",    "var":baseWeight},
{"name":"PU_Up",     "tree":"nom2",    "var":"(tri*weight*puweightUp*mueffweight*eleffweight)"},
#{"name":"PU_Down",   "tree":"nom2",    "var":"(tri*weight*puweightDown*mueffweight*eleffweight)"},
{"name":"JER_Up",    "tree":"nomJER_up",   "var":baseWeight},
#{"name":"JER_Down",  "tree":"nomJER_dw",   "var":baseWeight},
{"name":"JES_Up",    "tree":"nomJES_up",   "var":baseWeight},        
#{"name":"JES_Down",  "tree":"nomJES_dw",   "var":baseWeight},      
{"name":"MuEff_Up",     "tree":"nom2",    "var":"(tri*weight*puweight*mueffweight_up*eleffweight)"},
#{"name":"MuEff_Down",   "tree":"nom2",    "var":"(tri*weight*puweight*mueffweight_dn*eleffweight)"},
{"name":"ElEff_Up",     "tree":"nom2",    "var":"(tri*weight*puweight*mueffweight*eleffweight_up)"},
#{"name":"ElEff_Down",   "tree":"nom2",    "var":"(tri*weight*puweight*mueffweight*eleffweight_dn)"},
{"name":"Mu_Scale_Up",     "tree":"nomMu_up",    "var":"("+baseWeight+"*csvweights2[0])"},
#{"name":"Mu_Scale_Down",     "tree":"nomMu_dw",    "var":"("+baseWeight+"*csvweights2[0])"},
{"name":"El_Scale_Up",     "tree":"nomEl_up",    "var":"("+baseWeight+"*csvweights2[0])"},
#{"name":"El_Scale_Down",     "tree":"nomEl_dw",    "var":"("+baseWeight+"*csvweights2[0])"},

{"name":"Trigger_Up",   "tree":"nom2",        "var": "(weight*puweight*mueffweight*eleffweight*tri_up)"},
#{"name":"Trigger_Down", "tree":"nom2",        "var": "(weight*puweight*mueffweight*eleffweight*tri_dn)"},
{"name":"topPt_Up", "tree":"nom2",      "var":"("+baseWeight+"*topPtWeight)"},
#{"name":"topPt_Down", "tree":"nom2",    "var":"("+baseWeight+")"},
#{"name":"FSR_Down", "tree":"nom2",      "var":"("+baseWeight+")"},
{"name":"FSR_Up", "tree":"nom2",        "var":"("+baseWeight+")"},
#{"name":"ISR_Down", "tree":"nom2",      "var":"("+baseWeight+")"},
{"name":"ISR_Up", "tree":"nom2",        "var":"("+baseWeight+")"},
#{"name":"UE_Down", "tree":"nom2",       "var":"("+baseWeight+")"},
{"name":"UE_Up", "tree":"nom2",         "var":"("+baseWeight+")"},
#{"name":"hdamp_Down", "tree":"nom2",    "var":"("+baseWeight+")"},
{"name":"hdamp_Up", "tree":"nom2",      "var":"("+baseWeight+")"},
#{"name":"qcderdon_Down", "tree":"nom2", "var":"("+baseWeight+")"},
{"name":"qcderdon_Up", "tree":"nom2",   "var":"("+baseWeight+")"},
#{"name":"erdon_Down", "tree":"nom2",    "var":"("+baseWeight+")"},
{"name":"erdon_Up", "tree":"nom2",      "var":"("+baseWeight+")"},
#{"name":"gluonmove_Down", "tree":"nom2","var":"("+baseWeight+")"},
{"name":"gluonmove_Up", "tree":"nom2",  "var":"("+baseWeight+")"},
{"name":"MuF_Up", "tree":"nom2",        "var":"("+baseWeight+"*scaleWeightsUp[0])"},
#{"name":"MuF_Down", "tree":"nom2",      "var":"("+baseWeight+"*scaleWeightsDown[0])"},
#{"name":"MuR_Down", "tree":"nom2",      "var":"("+baseWeight+"*scaleWeightsDown[1])"},
{"name":"MuR_Up", "tree":"nom2",        "var":"("+baseWeight+"*scaleWeightsUp[1])"},
{"name":"PDF_Up", "tree":"nom2",        "var":"("+baseWeight+"*pdfWeights[0])"},
#{"name":"PDF_Down", "tree":"nom2",      "var":"("+baseWeight+"*pdfWeights[1])"},
{"name":"PDFAlphaS_Up", "tree":"nom2",  "var":"("+baseWeight+"*pdfWeights[2])"},
#{"name":"PDFAlphaS_Down", "tree":"nom2","var":"("+baseWeight+"*pdfWeights[3])"},
]
mceventweightCSV=[
#{"name":"CEN",       "tree":"nom2",    "var":"("+baseWeight+")"},
#{"name":"NoPU",      "tree":"nom2",    "var":"("+baseWeightNOPU+")"},
{"name":"csvweight", "tree":"nom2",    "var":baseWeight2},
{"name":"PU_Up",     "tree":"nom2",    "var":"(tri*weight*puweightUp*mueffweight*eleffweight*csvweights2[0])"},
#]
#test=[
#{"name":"PU_Down",   "tree":"nom2",    "var":"(tri*weight*puweightDown*mueffweight*eleffweight*csvweights2[0])"},
{"name":"JER_Up",    "tree":"nomJER_up",   "var":baseWeight2},
#{"name":"JER_Down",  "tree":"nomJER_dw",   "var":baseWeight2},
{"name":"JES_Up",    "tree":"nomJES_up",   "var":baseWeight2},        
#{"name":"JES_Down",  "tree":"nomJES_dw",   "var":baseWeight2},      
{"name":"MuEff_Up",     "tree":"nom2",    "var":"(tri*weight*puweight*mueffweight_up*eleffweight*csvweights2[0])"},
#{"name":"MuEff_Down",   "tree":"nom2",    "var":"(tri*weight*puweight*mueffweight_dn*eleffweight*csvweights2[0])"},
{"name":"ElEff_Up",     "tree":"nom2",    "var":"(tri*weight*puweight*mueffweight*eleffweight_up*csvweights2[0])"},
#{"name":"ElEff_Down",   "tree":"nom2",    "var":"(tri*weight*puweight*mueffweight*eleffweight_dn*csvweights2[0])"},
{"name":"Mu_Scale_Up",     "tree":"nomMu_up",    "var":"("+baseWeight+"*csvweights2[0])"},
#{"name":"Mu_Scale_Down",     "tree":"nomMu_dw",    "var":"("+baseWeight+"*csvweights2[0])"},
{"name":"El_Scale_Up",     "tree":"nomEl_up",    "var":"("+baseWeight+"*csvweights2[0])"},
#{"name":"El_Scale_Down",     "tree":"nomEl_dw",    "var":"("+baseWeight+"*csvweights2[0])"},

{"name":"Btag_LF_Up",     "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[3])"},
#{"name":"Btag_LF_Down",   "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[4])"},
{"name":"Btag_HF_Up",           "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[5])"}, 
#{"name":"Btag_HF_Down",         "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[6])"},
{"name":"Btag_HF_Stats1_Up",    "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[7])"},
#{"name":"Btag_HF_Stats1_Down",  "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[8])"},
{"name":"Btag_HF_Stats2_Up",    "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[9])"},
#{"name":"Btag_HF_Stats2_Down",  "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[10])"},
{"name":"Btag_LF_Stats1_Up",    "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[11])"},
#{"name":"Btag_LF_Stats1_Down",  "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[12])"},
{"name":"Btag_LF_Stats2_Up",    "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[13])"},
#{"name":"Btag_LF_Stats2_Down",  "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[14])"},
{"name":"Btag_CQ_Err1_Up",      "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[15])"},
#{"name":"Btag_CQ_Err1_Down",    "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[16])"},
{"name":"Btag_CQ_Err2_Up",      "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[17])"},
#{"name":"Btag_CQ_Err2_Down",    "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[18])"},
{"name":"Trigger_Up",   "tree":"nom2",        "var": "(weight*puweight*mueffweight*eleffweight*tri_up*csvweights2[0])"},
#{"name":"Trigger_Down", "tree":"nom2",        "var": "(weight*puweight*mueffweight*eleffweight*tri_dn*csvweights2[0])"},
{"name":"topPt_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0]*topPtWeight)"},
#{"name":"topPt_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
#{"name":"FSR_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
{"name":"FSR_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
#{"name":"ISR_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
{"name":"ISR_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
#{"name":"UE_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
{"name":"UE_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
#{"name":"hdamp_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
{"name":"hdamp_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
#{"name":"qcderdon_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
{"name":"qcderdon_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
#{"name":"erdon_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
{"name":"erdon_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
#{"name":"gluonmove_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
{"name":"gluonmove_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
{"name":"MuF_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0]*scaleWeightsUp[0])"},
#{"name":"MuF_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0]*scaleWeightsDown[0])"},
#{"name":"MuR_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0]*scaleWeightsDown[1])"},
{"name":"MuR_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0]*scaleWeightsUp[1])"},
{"name":"PDF_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0]*pdfWeights[0])"},
#{"name":"PDF_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0]*pdfWeights[1])"},
{"name":"PDFAlphaS_Up", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0]*pdfWeights[2])"},
#{"name":"PDFAlphaS_Down", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0]*pdfWeights[3])"},
]
#{"name":"Free_ttb_ttcc", "tree":"nom2",    "var":"("+baseWeight+"*csvweights2[0])"},
#mceventweightMG5 = [{'name':i['name'],'var':i['var'].replace("(weight*","(").replace('*weight*','*')} for i in mceventweight]
