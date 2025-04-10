# 011 "Scenarios"
# selection of input station data for the reconstruction and cross-validation. Scenario refers to the network that is chosen. 
# If a station does not exist, it is simply not used. 

all_temp = c("MIL_ta","TOR_ta","BOL_ta","BAS_ta","SMA_ta","BER_ta","HOH_ta","PAD_ta","SHA_ta","GVE_ta","STG_ta","ROV_ta","NEU_ta","FON_ta","BUS_ta","ZUG_ta","VEV_ta","GOT_ta",
             "ALT_ta","ANT_ta","RAG_ta","CHM_ta","CDF_ta","GSB_ta","DAV_ta","ELM_ta","ENG_ta","GRC_ta","MER_ta","OTL_ta","LUG_ta","LUZ_ta",
             "SAE_ta","SAM_ta","SBE_ta","SIA_ta","SIO_ta","CHD_ta","GRH_ta", "JUN_ta")

scen_all <- c("SBE_rr", "SVG_rr", "SOG_rr", "STB_rr", "SUS_rr", "THU_rr", "TST_rr", "WEE_rr", "AIE_rr", "AIR_rr", "ALS_rr", "ANT_rr", "APP_rr", "BEX_rr", 
              "BIN_rr", "BIV_rr", "BRP_rr", "CHD_rr", "GSB_rr", "COU_rr", "ELM_rr", "ENG_rr", "ENT_rr", "ESZ_rr", "FLU_rr", "GTT_rr", "GOS_rr", "ILZ_rr", 
              "KAS_rr", "KLA_rr", "CDF_rr", "VST_rr", "LAC_rr", "LAB_rr", "LSN_rr", "LTB_rr", "LEU_rr", "OTL_rr", "LOH_rr", "LON_rr", "MAT_rr", "MER_rr", 
              "MMO_rr", "MSG_rr", "MUR_rr", "OBI_rr", "ROM_rr", "SAR_rr", "SAE_rr", "BAS_rr", "BER_rr", "BIA_rr", "CHM_rr", "DAV_rr", "GVE_rr", "GRC_rr", 
              "LUG_rr", "LUZ_rr", "NEU_rr", "SAM_rr", "SIA_rr", "SIO_rr", "STG_rr", "WIN_rr", "WIT_rr", "SMA_rr", "ALT_rr", "THS_rr", "ALT_ta", "ANT_ta", 
              "RAG_ta", "BAS_ta", "BER_ta", "CHM_ta", "CHD_ta", "GSB_ta", "DAV_ta", "ELM_ta", "ENG_ta", "GVE_ta", "GRH_ta", "GRC_ta", "JUN_ta", "CDF_ta", 
              "OTL_ta", "LUG_ta", "LUZ_ta", "MER_ta", "NEU_ta", "SBE_ta", "SAM_ta", "SIA_ta", "SIO_ta", "STG_ta", "SAE_ta", "SMA_ta", "ALT_p", "BAS_p", 
              "BER_p", "CHU_p", "GSB_p", "DAV_p", "CDF_p", "LUG_p", "NEU_p", "SHA_p", "SAE_p", "SMA_p", "BAS_rr0", "BER_rr0", "BOL_ta", "BUS_p", "BUS_rr0", 
              "BUS_ta", "BUS_rr", "ENS_rr0", "FON_ta", "GOT_p", "GOT_ta", "GVE_p", "HOH_p", "HOH_ta", "KAR_p", "LUZ_p", "MIL_p", "MIL_ta", "PAD_p", "PAD_ta",
              "ROV_ta", "SHA_ta", "STG_rr0", "STG_p", "TOR_ta", "TOR_p", "VEV_ta", "ZUG_ta", "NSG_sp", "EWG_sp", "LUG_sp", "ALT_sp", "BAS_sp", "BER_sp", 
              "CHU_sp", "GSB_sp", "DAV_sp", "CDF_sp", "NEU_sp", "SHA_sp", "SAE_sp", "SMA_sp", "BUS_sp", "GOT_sp", "GVE_sp", "HOH_sp", "KAR_sp", "LUZ_sp", 
              "MIL_sp", "PAD_sp", "STG_sp", "TOR_sp", "MIL_spm", "TOR_spm", "ALT_spm", "BAS_spm", "BER_spm", "CHU_spm", "GSB_spm", "DAV_spm", 
              "CDF_spm", "LUG_spm", "NEU_spm", "SHA_spm", "SAE_spm", "SMA_spm", "BUS_spm", "GOT_spm", "GVE_spm", "HOH_spm", "KAR_spm", "LUZ_spm", "STG_spm")


scen_hist_wind <- c("MIL_p","BAS_p","SMA_p","NEU_p","BER_p","TOR_p","HOH_p","PAD_p","GOT_p","GVE_p","STG_p","BUS_p","KAR_p","GSB_p","SHA_p","LUZ_p",
                    "BAS_rr0","BER_rr","BUS_rr0","GVE_rr","STG_rr0","ENS_rr0","BUS_rr",
                    "MIL_ta","TOR_ta","BOL_ta","BAS_ta","SMA_ta","NEU_ta","BER_ta","HOH_ta","PAD_ta","GOT_ta",
                    "SHA_ta","GVE_ta","STG_ta","ROV_ta","FON_ta","BUS_ta","ZUG_ta","GSB_ta","VEV_ta","LUZ_ta",
                    "MIL_spm", "TOR_spm","BAS_spm", "BER_spm", "GSB_spm",
                    "NEU_spm", "SHA_spm", "SMA_spm", "BUS_spm", "GOT_spm", "GVE_spm", "HOH_spm", "KAR_spm", "LUZ_spm", "STG_spm",
                    "EWG_sp")

scen_hist2_wind = c("ALT_p","BAS_p","BER_p","CHU_p","GSB_p","DAV_p","CDF_p","LUG_p","NEU_p","SHA_p","SAE_p","SMA_p",
                    "ALT_ta","RAG_ta","BAS_ta","BER_ta","CHM_ta","GSB_ta","DAV_ta","ELM_ta","GVE_ta","GRH_ta","GRC_ta","JUN_ta",
                    "CDF_ta","OTL_ta","LUG_ta","MER_ta","SIA_ta","SIO_ta","SAE_ta","SMA_ta",
                    "AIR_rr0","SAE_rr0","GSB_rr0","LEU_rr0","CDF_rr0","LON_rr0","STB_rr0",
                    "BAS_spm","BER_spm","NEU_spm","SHA_spm","SMA_spm",
                    "EWG_sp") 

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


scen_sel = c("ALT_p","BAS_p","BER_p","CHU_p","GSB_p","DAV_p","CDF_p","LUG_p","NEU_p","SHA_p","SAE_p","SMA_p",
             "ALT_ta","RAG_ta","BAS_ta","BER_ta","CHM_ta","GSB_ta","DAV_ta","ELM_ta","GVE_ta","GRH_ta","GRC_ta","JUN_ta",
             "CDF_ta","OTL_ta","LUG_ta","MER_ta","SIA_ta","SIO_ta","SAE_ta","SMA_ta",
             "AIR_rr","SAE_rr","GSB_rr","LEU_rr","CDF_rr","LON_rr","STB_rr")


scen4 = c("ALT_ta", "ANT_ta", "RAG_ta", "BAS_ta", "BER_ta", "CHM_ta", "GSB_ta","GRC_ta" , "DAV_ta", "ELM_ta", "ENG_ta", "GVE_ta", 
          "LUZ_ta", "MER_ta" ,"NEU_ta", "SBE_ta",
          "SAM_ta", "SIA_ta", "SIO_ta", "STG_ta", "SAE_ta", "SMA_ta", "OTL_ta", "LUG_ta",
          "ALT_p",  "BAS_p",  "BER_p",  "CHU_p",  "GSB_p",  "DAV_p",  "LUG_p",  "NEU_p",  "SHA_p",  "SAE_p", "SMA_p",
          "SVG_rr", "SOG_rr", "THU_rr", "WEE_rr", "AIE_rr", "AIR_rr", "ALS_rr", "ANT_rr", "BIV_rr", "ELM_rr", "ENT_rr", "ESZ_rr", "FLU_rr", "GTT_rr", "ILZ_rr" ,"VST_rr", "LAC_rr", "LAB_rr" ,"LSN_rr","LEU_rr",
          "OTL_rr", "LOH_rr", "LON_rr", "MER_rr", "MMO_rr", "MUR_rr" ,"OBI_rr" ,"SAR_rr" ,"SAE_rr" ,"BAS_rr" ,"BER_rr" ,"CHM_rr" ,"DAV_rr" ,"GVE_rr" ,"LUG_rr" ,"LUZ_rr", "NEU_rr", "SAM_rr", "SIA_rr" ,"SIO_rr",
          "STG_rr", "WIN_rr" ,"WIT_rr", "SMA_rr", "ALT_rr" ,"THS_rr" )


scen3 = c("BAS_ta","BER_ta","GSB_ta","GVE_ta","STG_ta","SMA_ta","BOL_ta", "BUS_ta","HOH_ta","MIL_ta","PAD_ta","ROV_ta","SHA_ta" , "TOR_ta" ,"VEV_ta",
          "BAS_p","BER_p","GSB_p","SHA_p","SMA_p", "BUS_p","GVE_p","HOH_p", "KAR_p", "MIL_p" ,"PAD_p",  "STG_p" , "TOR_p" ,
          "GVE_rr","BUS_rr0","ENS_rr0","STG_rr0",
          "BAS_spm","BER_spm","GSB_spm","SHA_spm","SMA_spm", "BUS_spm","GVE_spm","STG_spm" ,"HOH_spm", "KAR_spm", "MIL_spm" , "TOR_spm", 
          "EWG_sp")


scen2 = c("MIL_ta","TOR_ta","BOL_ta","BAS_ta","SMA_ta","BER_ta","HOH_ta","PAD_ta","GVE_ta","ROV_ta","SHA_ta",
          "MIL_p","BAS_p","BER_p","HOH_p","PAD_p","SHA_p","KAR_p","TOR_p",
          "GVE_rr","BAS_rr0",
          "MIL_spm","BAS_spm","BER_spm","HOH_spm","SHA_spm","KAR_spm","TOR_spm")

scen1 = c("MIL_ta","TOR_ta","BOL_ta","BAS_ta","SMA_ta","NEU_ta","PAD_ta",
          "MIL_p","BAS_p","SMA_p","NEU_p",
          "BAS_rr0",
          "MIL_spm","BAS_spm","SMA_spm","NEU_spm")



scen_spsel_rr0 = c("ALT_p","BAS_p","BER_p","CHU_p","GSB_p","DAV_p","CDF_p","LUG_p","NEU_p","SHA_p","SAE_p","SMA_p",
                   "ALT_ta","RAG_ta","BAS_ta","BER_ta","CHM_ta","GSB_ta","DAV_ta","ELM_ta","GVE_ta","GRH_ta","GRC_ta","JUN_ta",
                   "CDF_ta","OTL_ta","LUG_ta","MER_ta","SIA_ta","SIO_ta","SAE_ta","SMA_ta",
                   "AIR_rr0","SAE_rr0","GSB_rr0","LEU_rr0","CDF_rr0","LON_rr0","STB_rr0",
                   "BAS_sp","BER_sp","CHU_sp","GSB_sp","DAV_sp","CDF_sp","LUG_sp","NEU_sp","SHA_sp","SAE_sp","SMA_sp",
                   "EWG_sp")


scenarios <- list(scen_all_ta_p_rr = scen_all_ta_p_rr,
                  scen_all = scen_all,
                  scen4 = scen4,
                  scen3 = scen3,
                  scen2 = scen2,
                  scen1 = scen1,
                  scen_spsel_rr0 = scen_spsel_rr0,
                  scen_hist_wind = scen_hist_wind,
                  scen_hist2_wind = scen_hist2_wind,
                  all_temp = all_temp) 


