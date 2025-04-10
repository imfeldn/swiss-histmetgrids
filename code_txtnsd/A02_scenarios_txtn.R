#  "Scenarios" of station availability
# selection of input station data for the reconstruction and cross-validation. Scenario refers to the network that is chosen. 
# if a station does not exist, it is simply not used. 


scen5_tmax2 <- c( "ALT_ta" ,"ANT_ta", "RAG_ta", "BAS_ta", "BER_ta", "CHM_ta", "CHD_ta" ,"GSB_ta", "DAV_ta", "ELM_ta", "ENG_ta", "GVE_ta",
                  "GRH_ta", "GRC_ta", "JUN_ta", "CDF_ta", "OTL_ta", "LUG_ta", "LUZ_ta", "MER_ta", "NEU_ta", "SBE_ta", "SAM_ta", "SIA_ta",
                  "SIO_ta", "STG_ta", "SAE_ta", "SMA_ta", "ALT_p",  "BAS_p",  "BER_p",  "CHU_p" ,
                  "GSB_p",  "DAV_p",  "CDF_p",  "LUG_p" , "NEU_p" , "SHA_p" , "SAE_p" , "SMA_p",
                  "ALT_tmin","BAS_tmax" ,"BAS_tmin","BER_tmax","BER_tmin","DAV_tmax","DAV_tmin", "GVE_tmax","GVE_tmin","CDF_tmin", "OTL_tmin", 
                  "LUG_tmax", "LUG_tmin" ,"LUZ_tmax", "LUZ_tmin",     
                  "NEU_tmin","NEU_tmax", "SAM_tmax","SAM_tmin", "STG_tmin","SAE_tmax","SAE_tmin" ,"SMA_tmax","SMA_tmin")

scen5_tmax1 = c( "ALT_ta" ,"ANT_ta", "RAG_ta", "BAS_ta", "BER_ta", "CHM_ta", "CHD_ta" ,"GSB_ta", "DAV_ta", "ELM_ta", "ENG_ta", "GVE_ta",
                 "GRH_ta", "GRC_ta", "JUN_ta", "CDF_ta", "OTL_ta", "LUG_ta", "LUZ_ta", "MER_ta", "NEU_ta", "SBE_ta", "SAM_ta", "SIA_ta",
                 "SIO_ta", "STG_ta", "SAE_ta", "SMA_ta", "ALT_p",  "BAS_p",  "BER_p",  "CHU_p" ,
                 "GSB_p",  "DAV_p",  "CDF_p",  "LUG_p" , "NEU_p" , "SHA_p" , "SAE_p" , "SMA_p" )


scen_all_ta_p_rr <- c("SBE_rr", "SVG_rr", "SOG_rr", "STB_rr", "SUS_rr", "THU_rr", "TST_rr", "WEE_rr", "AIE_rr", "AIR_rr", "ALS_rr", "ANT_rr", "APP_rr", "BEX_rr", 
                      "BIN_rr", "BIV_rr", "BRP_rr", "CHD_rr", "GSB_rr", "COU_rr", "ELM_rr", "ENG_rr", "ENT_rr", "ESZ_rr", "FLU_rr", "GTT_rr", "GOS_rr", "ILZ_rr", 
                      "KAS_rr", "KLA_rr", "CDF_rr", "VST_rr", "LAC_rr", "LAB_rr", "LSN_rr", "LTB_rr", "LEU_rr", "OTL_rr", "LOH_rr", "LON_rr", "MAT_rr", "MER_rr", 
                      "MMO_rr", "MSG_rr", "MUR_rr", "OBI_rr", "ROM_rr", "SAR_rr", "SAE_rr", "BAS_rr", "BER_rr", "BIA_rr", "CHM_rr", "DAV_rr", "GVE_rr", "GRC_rr", 
                      "LUG_rr", "LUZ_rr", "NEU_rr", "SAM_rr", "SIA_rr", "SIO_rr", "STG_rr", "WIN_rr", "WIT_rr", "SMA_rr", "ALT_rr", "THS_rr", "ALT_ta", "ANT_ta", 
                      "RAG_ta", "BAS_ta", "BER_ta", "CHM_ta", "CHD_ta", "GSB_ta", "DAV_ta", "ELM_ta", "ENG_ta", "GVE_ta", "GRH_ta", "GRC_ta", "JUN_ta", "CDF_ta", 
                      "OTL_ta", "LUG_ta", "LUZ_ta", "MER_ta", "NEU_ta", "SBE_ta", "SAM_ta", "SIA_ta", "SIO_ta", "STG_ta", "SAE_ta", "SMA_ta", "ALT_p", "BAS_p", 
                      "BER_p", "CHU_p", "GSB_p", "DAV_p", "CDF_p", "LUG_p", "NEU_p", "SHA_p", "SAE_p", "SMA_p", "BAS_rr0", "BER_rr0", "BOL_ta", "BUS_p", "BUS_rr0", 
                      "BUS_ta", "BUS_rr", "ENS_rr0", "FON_ta", "GOT_p", "GOT_ta", "GVE_p", "HOH_p", "HOH_ta", "KAR_p", "LUZ_p", "MIL_p", "MIL_ta", "PAD_p", "PAD_ta",
                      "ROV_ta", "SHA_ta", "STG_rr0", "STG_p", "TOR_ta", "TOR_p", "VEV_ta", "ZUG_ta")


scen4_tmax2 <- c( "ALT_ta" ,"ANT_ta", "RAG_ta", "BAS_ta", "BER_ta", "CHM_ta", "GSB_ta", "GRC_ta","DAV_ta", "ELM_ta", "ENG_ta", "GVE_ta",
                  "OTL_ta", "LUG_ta", "LUZ_ta", "MER_ta", "NEU_ta", "SBE_ta", "SAM_ta", "SIA_ta",
                  "STG_ta", "SAE_ta", "SMA_ta",
                  "ALT_p",  "BAS_p",  "BER_p",  "CHU_p" ,"GSB_p",  "DAV_p",  "LUG_p" , "NEU_p" , "SHA_p" , "SAE_p" , "SMA_p",
                  "BER_tmax","BER_tmin","GVE_tmax","GVE_tmin",
                  "LUG_tmax", "LUG_tmin","NEU_tmin","NEU_tmax")

scen3 = c("BAS_ta","BER_ta","GSB_ta","GVE_ta","STG_ta","SMA_ta","BOL_ta", "BUS_ta","HOH_ta","MIL_ta","PAD_ta","ROV_ta","SHA_ta" , "TOR_ta" ,"VEV_ta",
          "BAS_p","BER_p","GSB_p","SHA_p","SMA_p", "BUS_p","GVE_p","HOH_p", "KAR_p", "MIL_p" ,"PAD_p",  "STG_p" , "TOR_p" ,
          "GVE_rr","BUS_rr0","ENS_rr0","STG_rr0")

scen1 = c("MIL_ta","TOR_ta","BOL_ta","BAS_ta","SMA_ta","NEU_ta",
          "MIL_p","BAS_p","SMA_p","NEU_p","PAD_ta",
          "BAS_rr0")

scen2 = c("MIL_ta","TOR_ta","BOL_ta","BAS_ta","SMA_ta","BER_ta","HOH_ta","PAD_ta","GVE_ta","ROV_ta","SHA_ta",
          "MIL_p","BAS_p","BER_p","HOH_p","PAD_p","SHA_p","KAR_p","TOR_p",
          "GVE_rr","BAS_rr0")


all_hist = unique(c(scen1,scen2,scen3, "FON_ta","ZUG_ta","GOT_ta")) # add the stations in none of the selected networks examples
all_hist2 = unique(c(scen4_tmax2,scen5_tmax2)) # add the stations in none of the selected networks examples
all_temp_txtn = all_hist[grepl("_ta|_tmax|_tmin",all_hist)]

scenarios <- list(all_temp_txtn = all_temp_txtn, 
                  all_hist = all_hist,
                  all_hist2 = all_hist2,
                  scen_all_ta_p_rr = scen_all_ta_p_rr,
                  scen5_tmax1 = scen5_tmax1,
                  scen5_tmax2 = scen5_tmax2,
                  scen4_tmax2 = scen4_tmax2,
                  scen3 = scen3,
                  scen2 = scen2,
                  scen1= scen1
) 

# 