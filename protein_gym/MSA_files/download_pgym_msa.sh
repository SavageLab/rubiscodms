#!/bin/bash

# Check if a path is provided
if [ "$#" -ne 1 ]; then
    echo "Please provide an output_path where MSAs will be stored: $0 <output_path>"
    exit 1
fi

# Assign the provided path to a variable
output_path="$1"

curl -o MSA_ProteinGym.zip https://marks.hms.harvard.edu/proteingym/MSA_ProteinGym.zip
unzip MSA_ProteinGym.zip
rm MSA_ProteinGym.zip

# Define a function to copy files to the specified directory
copy_and_create_directory() {
    local source_file="$1"
    local target_directory="$2"
    mkdir -p "$output_path/$target_directory"
    cp "$source_file" "$output_path/$target_directory/nt_alignment_msa.fna"
}

# Copy files to the specified directories
copy_and_create_directory MSA_files/A0A140D2T1_ZIKV_theta0.99_281-804_11-26-2021_b02.a2m A0A140D2T1_ZIKV_Sourisseau_growth_2019/processed_inputs
copy_and_create_directory MSA_files/A0A192B1T2_9HIV1_theta0.99_full_11-26-2021_b09.a2m A0A192B1T2_9HIV1_Haddox_2018/processed_inputs
copy_and_create_directory MSA_files/A0A1I9GEU1_NEIME_full_11-26-2021_b08.a2m A0A1I9GEU1_NEIME_Kennouche_2019/processed_inputs
copy_and_create_directory MSA_files/A0A2Z5U3Z0_9INFA_theta0.99_full_11-26-2021_b09.a2m A0A2Z5U3Z0_9INFA_Doud_2016/processed_inputs
copy_and_create_directory MSA_files/A0A2Z5U3Z0_9INFA_theta0.99_full_11-26-2021_b09.a2m A0A2Z5U3Z0_9INFA_Wu_2014/processed_inputs
copy_and_create_directory MSA_files/A4_HUMAN_full_11-26-2021_b01.a2m A4_HUMAN_Seuma_2021/processed_inputs
copy_and_create_directory MSA_files/A4D664_9INFA_theta0.99_full_11-26-2021_b09.a2m A4D664_9INFA_Soh_CCL141_2019/processed_inputs
copy_and_create_directory MSA_files/A4GRB6_PSEAI_full_11-26-2021_b03.a2m A4GRB6_PSEAI_Chen_2020/processed_inputs
copy_and_create_directory MSA_files/AACC1_PSEAI_full_04-29-2022_b03.a2m AACC1_PSEAI_Dandage_2018/processed_inputs
copy_and_create_directory MSA_files/ADRB2_HUMAN_full_11-26-2021_b03.a2m ADRB2_HUMAN_Jones_2020/processed_inputs
copy_and_create_directory MSA_files/AMIE_PSEAE_full_11-26-2021_b02.a2m AMIE_PSEAE_Wrenbeck_2017/processed_inputs
copy_and_create_directory MSA_files/B3VI55_LIPST_full_11-26-2021_b03.a2m B3VI55_LIPST_Klesmith_2015/processed_inputs
copy_and_create_directory MSA_files/BLAT_ECOLX_full_11-26-2021_b02.a2m BLAT_ECOLX_Deng_2012/processed_inputs
copy_and_create_directory MSA_files/BLAT_ECOLX_full_11-26-2021_b02.a2m BLAT_ECOLX_Firnberg_2014/processed_inputs
copy_and_create_directory MSA_files/BLAT_ECOLX_full_11-26-2021_b02.a2m BLAT_ECOLX_Jacquier_2013/processed_inputs
copy_and_create_directory MSA_files/BLAT_ECOLX_full_11-26-2021_b02.a2m BLAT_ECOLX_Stiffler_2015/processed_inputs
copy_and_create_directory MSA_files/BRCA1_HUMAN_full_11-26-2021_b02.a2m BRCA1_HUMAN_Findlay_2018/processed_inputs
copy_and_create_directory MSA_files/C6KNH7_9INFA_theta0.99_full_11-26-2021_b09.a2m C6KNH7_9INFA_Lee_2018/processed_inputs
copy_and_create_directory MSA_files/CALM1_HUMAN_full_11-26-2021_b03.a2m CALM1_HUMAN_Weile_2017/processed_inputs
copy_and_create_directory MSA_files/CAPSD_AAV2S_uniprot_t099_msc70_mcc70_b0.8.a2m CAPSD_AAV2S_Sinai_substitutions_2021/processed_inputs
copy_and_create_directory MSA_files/CCDB_ECOLI_full_11-26-2021_b02.a2m CCDB_ECOLI_Adkar_2012/processed_inputs
copy_and_create_directory MSA_files/CCDB_ECOLI_full_11-26-2021_b02.a2m CCDB_ECOLI_Tripathi_2016/processed_inputs
copy_and_create_directory MSA_files/CP2C9_HUMAN_full_11-26-2021_b04.a2m CP2C9_HUMAN_Amorosi_abundance_2021/processed_inputs
copy_and_create_directory MSA_files/CP2C9_HUMAN_full_11-26-2021_b04.a2m CP2C9_HUMAN_Amorosi_activity_2021/processed_inputs
copy_and_create_directory MSA_files/DLG4_HUMAN_full_11-26-2021_b02.a2m DLG4_HUMAN_Faure_2021/processed_inputs
copy_and_create_directory MSA_files/DLG4_RAT_full_11-26-2021_b03.a2m DLG4_RAT_McLaughlin_2012/processed_inputs
copy_and_create_directory MSA_files/DYR_ECOLI_full_11-26-2021_b08.a2m DYR_ECOLI_Thompson_plusLon_2019/processed_inputs
copy_and_create_directory MSA_files/ENV_HV1B9_S364P-M373R_b0.3.a2m ENV_HV1B9_DuenasDecamp_2016/processed_inputs
copy_and_create_directory MSA_files/ENV_HV1BR_theta0.99_full_11-26-2021_b09.a2m ENV_HV1BR_Haddox_2016/processed_inputs
copy_and_create_directory MSA_files/ESTA_BACSU_full_11-26-2021_b03.a2m ESTA_BACSU_Nutschel_2020/processed_inputs
copy_and_create_directory MSA_files/F7YBW8_MESOW_full_01-07-2022_b02.a2m F7YBW8_MESOW_Aakre_2015/processed_inputs
copy_and_create_directory MSA_files/GAL4_YEAST_full_11-26-2021_b02.a2m GAL4_YEAST_Kitzman_2015/processed_inputs
copy_and_create_directory MSA_files/GCN4_YEAST_full_24-02-2022_b03.a2m GCN4_YEAST_Staller_induction_2018/processed_inputs
copy_and_create_directory MSA_files/GFP_AEQVI_full_04-29-2022_b08.a2m GFP_AEQVI_Sarkisyan_2016/processed_inputs
copy_and_create_directory MSA_files/GRB2_HUMAN_full_11-26-2021_b05.a2m GRB2_HUMAN_Faure_2021/processed_inputs
copy_and_create_directory MSA_files/HIS7_YEAST_full_11-26-2021_b09.a2m HIS7_YEAST_Pokusaeva_2019/processed_inputs
copy_and_create_directory MSA_files/HSP82_YEAST_full_11-26-2021_b01.a2m HSP82_YEAST_Flynn_2019/processed_inputs
copy_and_create_directory MSA_files/HSP82_YEAST_full_11-26-2021_b01.a2m HSP82_YEAST_Mishra_2016/processed_inputs
copy_and_create_directory MSA_files/I6TAH8_I68A0_theta0.99_full_11-26-2021_b09.a2m I6TAH8_I68A0_Doud_2015/processed_inputs
copy_and_create_directory MSA_files/IF1_ECOLI_full_11-26-2021_b02.a2m IF1_ECOLI_Kelsic_2016/processed_inputs
copy_and_create_directory MSA_files/KCNH2_HUMAN_535-565_11-26-2021_b05.a2m KCNH2_HUMAN_Kozek_2020/processed_inputs
copy_and_create_directory MSA_files/KKA2_KLEPN_full_11-26-2021_b02.a2m KKA2_KLEPN_Melnikov_2014/processed_inputs
copy_and_create_directory MSA_files/MK01_HUMAN_full_11-26-2021_b06.a2m MK01_HUMAN_Brenan_2016/processed_inputs
copy_and_create_directory MSA_files/MSH2_HUMAN_full_11-26-2021_b05.a2m MSH2_HUMAN_Jia_2020/processed_inputs
copy_and_create_directory MSA_files/MTH3_HAEAE_full_11-26-2021_b02.a2m MTH3_HAEAE_Rockah-Shmuel_2015/processed_inputs
copy_and_create_directory MSA_files/NCAP_I34A1_theta0.99_full_11-26-2021_b09.a2m NCAP_I34A1_Doud_2015/processed_inputs
copy_and_create_directory MSA_files/NRAM_I33A0_full_11-26-2021_b01.a2m NRAM_I33A0_Jiang_standard_2016/processed_inputs
copy_and_create_directory MSA_files/NUD15_HUMAN_full_11-26-2021_b04.a2m NUD15_HUMAN_Suiter_2020/processed_inputs
copy_and_create_directory MSA_files/P53_HUMAN_full_04-29-2022_b09.a2m P53_HUMAN_Giacomelli_NULL_Etoposide_2018/processed_inputs
copy_and_create_directory MSA_files/P53_HUMAN_full_04-29-2022_b09.a2m P53_HUMAN_Giacomelli_NULL_Nutlin_2018/processed_inputs
copy_and_create_directory MSA_files/P53_HUMAN_full_04-29-2022_b09.a2m P53_HUMAN_Giacomelli_WT_Nutlin_2018/processed_inputs
copy_and_create_directory MSA_files/P53_HUMAN_full_11-26-2021_b09.a2m P53_HUMAN_Kotler_2018/processed_inputs
copy_and_create_directory MSA_files/P84126_THETH_full_11-26-2021_b04.a2m P84126_THETH_Chan_2017/processed_inputs
copy_and_create_directory MSA_files/PA_I34A1_full_theta0.99_04-29-2022_b09.a2m PA_I34A1_Wu_2015/processed_inputs
copy_and_create_directory MSA_files/PABP_YEAST_full_11-26-2021_b07.a2m PABP_YEAST_Melamed_2013/processed_inputs
copy_and_create_directory MSA_files/POLG_CXB3N_1-861_theta0.99_04-29-2022_b07.a2m POLG_CXB3N_Mattenberger_2021/processed_inputs
copy_and_create_directory MSA_files/POLG_HCVJF_theta0.99_1984-2089_11-26-2021_b08.a2m POLG_HCVJF_Qi_2014/processed_inputs
copy_and_create_directory MSA_files/PTEN_HUMAN_full_11-26-2021_b01.a2m PTEN_HUMAN_Matreyek_2021/processed_inputs
copy_and_create_directory MSA_files/PTEN_HUMAN_full_11-26-2021_b01.a2m PTEN_HUMAN_Mighell_2018/processed_inputs
copy_and_create_directory MSA_files/Q2N0S5_9HIV1_full_theta0.99_04-29-2022_b09.a2m Q2N0S5_9HIV1_Haddox_2018/processed_inputs
copy_and_create_directory MSA_files/Q59976_STRSQ_full_11-26-2021_b03.a2m Q59976_STRSQ_Romero_2015/processed_inputs
copy_and_create_directory MSA_files/R1AB_SARS2_02-19-2022_b07.a2m R1AB_SARS2_Flynn_growth_2022/processed_inputs
copy_and_create_directory MSA_files/RASH_HUMAN_full_11-26-2021_b03.a2m RASH_HUMAN_Bandaru_2017/processed_inputs
copy_and_create_directory MSA_files/REV_HV1H2_full_theta0.99_04-29-2022_b09.a2m REV_HV1H2_Fernandes_2016/processed_inputs
copy_and_create_directory MSA_files/RL401_YEAST_full_11-26-2021_b01.a2m RL401_YEAST_Mavor_2016/processed_inputs
copy_and_create_directory MSA_files/RL401_YEAST_full_11-26-2021_b01.a2m RL401_YEAST_Roscoe_2013/processed_inputs
copy_and_create_directory MSA_files/RL401_YEAST_full_11-26-2021_b01.a2m RL401_YEAST_Roscoe_2014/processed_inputs
copy_and_create_directory MSA_files/SC6A4_HUMAN_full_11-26-2021_b02.a2m SC6A4_HUMAN_Young_2021/processed_inputs
copy_and_create_directory MSA_files/SCN5A_HUMAN_1611-1642_11-26-2021_b03.a2m SCN5A_HUMAN_Glazer_2019/processed_inputs
copy_and_create_directory MSA_files/SPG1_STRSG_full_11-26-2021_b07.a2m SPG1_STRSG_Olson_2014/processed_inputs
copy_and_create_directory MSA_files/SPIKE_SARS2_theta0.99_full_11-26-2021_b01.a2m SPIKE_SARS2_Starr_bind_2020/processed_inputs
copy_and_create_directory MSA_files/SPIKE_SARS2_theta0.99_full_11-26-2021_b01.a2m SPIKE_SARS2_Starr_expr_2020/processed_inputs
copy_and_create_directory MSA_files/SRC_HUMAN_full_11-26-2021_b06.a2m SRC_HUMAN_Ahler_CD_2019/processed_inputs
copy_and_create_directory MSA_files/SUMO1_HUMAN_full_11-26-2021_b02.a2m SUMO1_HUMAN_Weile_2017/processed_inputs
copy_and_create_directory MSA_files/SYUA_HUMAN_full_04-29-2022_b01.a2m SYUA_HUMAN_Newberry_2020/processed_inputs
copy_and_create_directory MSA_files/TADBP_HUMAN_full_11-26-2021_b09.a2m TADBP_HUMAN_Bolognesi_2019/processed_inputs
copy_and_create_directory MSA_files/TAT_HV1BR_full_theta0.99_04-29-2022_b09.a2m TAT_HV1BR_Fernandes_2016/processed_inputs
copy_and_create_directory MSA_files/TPK1_HUMAN_full_11-26-2021_b02.a2m TPK1_HUMAN_Weile_2017/processed_inputs
copy_and_create_directory MSA_files/TPMT_HUMAN_full_11-26-2021_b03.a2m TPMT_HUMAN_Matreyek_2018/processed_inputs
copy_and_create_directory MSA_files/TPOR_HUMAN_full_11-26-2021_b01.a2m TPOR_HUMAN_Bridgford_S505N_2020/processed_inputs
copy_and_create_directory MSA_files/TRPC_SACS2_full_11-26-2021_b07.a2m TRPC_SACS2_Chan_2017/processed_inputs
copy_and_create_directory MSA_files/TRPC_THEMA_full_11-26-2021_b07.a2m TRPC_THEMA_Chan_2017/processed_inputs
copy_and_create_directory MSA_files/UBC9_HUMAN_full_11-26-2021_b03.a2m UBC9_HUMAN_Weile_2017/processed_inputs
copy_and_create_directory MSA_files/UBE4B_MOUSE_full_11-26-2021_b05.a2m UBE4B_MOUSE_Starita_2013/processed_inputs
copy_and_create_directory MSA_files/VKOR1_HUMAN_full_11-26-2021_b03.a2m VKOR1_HUMAN_Chiasson_abundance_2020/processed_inputs
copy_and_create_directory MSA_files/VKOR1_HUMAN_full_11-26-2021_b03.a2m VKOR1_HUMAN_Chiasson_activity_2020/processed_inputs
copy_and_create_directory MSA_files/YAP1_HUMAN_full_11-26-2021_b02.a2m YAP1_HUMAN_Araya_2012/processed_inputs
#
#
#
#curl -o MSA_ProteinGym.zip https://marks.hms.harvard.edu/proteingym/MSA_ProteinGym.zip
#unzip MSA_ProteinGym.zip
#rm MSA_ProteinGym.zip
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A0A140D2T1_ZIKV_Sourisseau_growth_2019/processed_inputs/
#cp A0A140D2T1_ZIKV_theta0.99_281-804_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A0A140D2T1_ZIKV_Sourisseau_growth_2019/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A0A192B1T2_9HIV1_Haddox_2018/processed_inputs/
#cp A0A192B1T2_9HIV1_theta0.99_full_11-26-2021_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A0A192B1T2_9HIV1_Haddox_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A0A1I9GEU1_NEIME_Kennouche_2019/processed_inputs/
#cp A0A1I9GEU1_NEIME_full_11-26-2021_b08.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A0A1I9GEU1_NEIME_Kennouche_2019/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A0A2Z5U3Z0_9INFA_Doud_2016/processed_inputs/
#cp A0A2Z5U3Z0_9INFA_theta0.99_full_11-26-2021_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A0A2Z5U3Z0_9INFA_Doud_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A0A2Z5U3Z0_9INFA_Wu_2014/processed_inputs/
#cp A0A2Z5U3Z0_9INFA_theta0.99_full_11-26-2021_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A0A2Z5U3Z0_9INFA_Wu_2014/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A4_HUMAN_Seuma_2021/processed_inputs/
#cp A4_HUMAN_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A4_HUMAN_Seuma_2021/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A4D664_9INFA_Soh_CCL141_2019/processed_inputs/
#cp A4D664_9INFA_theta0.99_full_11-26-2021_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A4D664_9INFA_Soh_CCL141_2019/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A4GRB6_PSEAI_Chen_2020/processed_inputs/
#cp A4GRB6_PSEAI_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/A4GRB6_PSEAI_Chen_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/AACC1_PSEAI_Dandage_2018/processed_inputs/
#cp AACC1_PSEAI_full_04-29-2022_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/AACC1_PSEAI_Dandage_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/ADRB2_HUMAN_Jones_2020/processed_inputs/
#cp ADRB2_HUMAN_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/ADRB2_HUMAN_Jones_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/AMIE_PSEAE_Wrenbeck_2017/processed_inputs/
#cp AMIE_PSEAE_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/AMIE_PSEAE_Wrenbeck_2017/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/B3VI55_LIPST_Klesmith_2015/processed_inputs/
#cp B3VI55_LIPST_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/B3VI55_LIPST_Klesmith_2015/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/BLAT_ECOLX_Deng_2012/processed_inputs/
#cp BLAT_ECOLX_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/BLAT_ECOLX_Deng_2012/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/BLAT_ECOLX_Firnberg_2014/processed_inputs/
#cp BLAT_ECOLX_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/BLAT_ECOLX_Firnberg_2014/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/BLAT_ECOLX_Jacquier_2013/processed_inputs/
#cp BLAT_ECOLX_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/BLAT_ECOLX_Jacquier_2013/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/BLAT_ECOLX_Stiffler_2015/processed_inputs/
#cp BLAT_ECOLX_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/BLAT_ECOLX_Stiffler_2015/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/BRCA1_HUMAN_Findlay_2018/processed_inputs/
#cp BRCA1_HUMAN_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/BRCA1_HUMAN_Findlay_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/C6KNH7_9INFA_Lee_2018/processed_inputs/
#cp C6KNH7_9INFA_theta0.99_full_11-26-2021_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/C6KNH7_9INFA_Lee_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CALM1_HUMAN_Weile_2017/processed_inputs/
#cp CALM1_HUMAN_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CALM1_HUMAN_Weile_2017/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CAPSD_AAV2S_Sinai_substitutions_2021/processed_inputs/
#cp CAPSD_AAV2S_uniprot_t099_msc70_mcc70_b0.8.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CAPSD_AAV2S_Sinai_substitutions_2021/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CCDB_ECOLI_Adkar_2012/processed_inputs/
#cp CCDB_ECOLI_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CCDB_ECOLI_Adkar_2012/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CCDB_ECOLI_Tripathi_2016/processed_inputs/
#cp CCDB_ECOLI_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CCDB_ECOLI_Tripathi_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CP2C9_HUMAN_Amorosi_abundance_2021/processed_inputs/
#cp CP2C9_HUMAN_full_11-26-2021_b04.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CP2C9_HUMAN_Amorosi_abundance_2021/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CP2C9_HUMAN_Amorosi_activity_2021/processed_inputs/
#cp CP2C9_HUMAN_full_11-26-2021_b04.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/CP2C9_HUMAN_Amorosi_activity_2021/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/DLG4_HUMAN_Faure_2021/processed_inputs/
#cp DLG4_HUMAN_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/DLG4_HUMAN_Faure_2021/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/DLG4_RAT_McLaughlin_2012/processed_inputs/
#cp DLG4_RAT_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/DLG4_RAT_McLaughlin_2012/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/DYR_ECOLI_Thompson_plusLon_2019/processed_inputs/
#cp DYR_ECOLI_full_11-26-2021_b08.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/DYR_ECOLI_Thompson_plusLon_2019/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/ENV_HV1B9_DuenasDecamp_2016/processed_inputs/
#cp ENV_HV1B9_S364P-M373R_b0.3.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/ENV_HV1B9_DuenasDecamp_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/ENV_HV1BR_Haddox_2016/processed_inputs/
#cp ENV_HV1BR_theta0.99_full_11-26-2021_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/ENV_HV1BR_Haddox_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/ESTA_BACSU_Nutschel_2020/processed_inputs/
#cp ESTA_BACSU_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/ESTA_BACSU_Nutschel_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/F7YBW8_MESOW_Aakre_2015/processed_inputs/
#cp F7YBW8_MESOW_full_01-07-2022_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/F7YBW8_MESOW_Aakre_2015/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/GAL4_YEAST_Kitzman_2015/processed_inputs/
#cp GAL4_YEAST_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/GAL4_YEAST_Kitzman_2015/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/GCN4_YEAST_Staller_induction_2018/processed_inputs/
#cp GCN4_YEAST_full_24-02-2022_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/GCN4_YEAST_Staller_induction_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/GFP_AEQVI_Sarkisyan_2016/processed_inputs/
#cp GFP_AEQVI_full_04-29-2022_b08.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/GFP_AEQVI_Sarkisyan_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/GRB2_HUMAN_Faure_2021/processed_inputs/
#cp GRB2_HUMAN_full_11-26-2021_b05.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/GRB2_HUMAN_Faure_2021/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/HIS7_YEAST_Pokusaeva_2019/processed_inputs/
#cp HIS7_YEAST_full_11-26-2021_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/HIS7_YEAST_Pokusaeva_2019/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/HSP82_YEAST_Flynn_2019/processed_inputs/
#cp HSP82_YEAST_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/HSP82_YEAST_Flynn_2019/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/HSP82_YEAST_Mishra_2016/processed_inputs/
#cp HSP82_YEAST_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/HSP82_YEAST_Mishra_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/I6TAH8_I68A0_Doud_2015/processed_inputs/
#cp I6TAH8_I68A0_theta0.99_full_11-26-2021_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/I6TAH8_I68A0_Doud_2015/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/IF1_ECOLI_Kelsic_2016/processed_inputs/
#cp IF1_ECOLI_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/IF1_ECOLI_Kelsic_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/KCNH2_HUMAN_Kozek_2020/processed_inputs/
#cp KCNH2_HUMAN_535-565_11-26-2021_b05.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/KCNH2_HUMAN_Kozek_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/KKA2_KLEPN_Melnikov_2014/processed_inputs/
#cp KKA2_KLEPN_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/KKA2_KLEPN_Melnikov_2014/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/MK01_HUMAN_Brenan_2016/processed_inputs/
#cp MK01_HUMAN_full_11-26-2021_b06.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/MK01_HUMAN_Brenan_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/MSH2_HUMAN_Jia_2020/processed_inputs/
#cp MSH2_HUMAN_full_11-26-2021_b05.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/MSH2_HUMAN_Jia_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/MTH3_HAEAE_Rockah-Shmuel_2015/processed_inputs/
#cp MTH3_HAEAE_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/MTH3_HAEAE_Rockah-Shmuel_2015/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/NCAP_I34A1_Doud_2015/processed_inputs/
#cp NCAP_I34A1_theta0.99_full_11-26-2021_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/NCAP_I34A1_Doud_2015/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/NRAM_I33A0_Jiang_standard_2016/processed_inputs/
#cp NRAM_I33A0_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/NRAM_I33A0_Jiang_standard_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/NUD15_HUMAN_Suiter_2020/processed_inputs/
#cp NUD15_HUMAN_full_11-26-2021_b04.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/NUD15_HUMAN_Suiter_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/P53_HUMAN_Giacomelli_NULL_Etoposide_2018/processed_inputs/
#cp P53_HUMAN_full_04-29-2022_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/P53_HUMAN_Giacomelli_NULL_Etoposide_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/P53_HUMAN_Giacomelli_NULL_Nutlin_2018/processed_inputs/
#cp P53_HUMAN_full_04-29-2022_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/P53_HUMAN_Giacomelli_NULL_Nutlin_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/P53_HUMAN_Giacomelli_WT_Nutlin_2018/processed_inputs/
#cp P53_HUMAN_full_04-29-2022_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/P53_HUMAN_Giacomelli_WT_Nutlin_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/P53_HUMAN_Kotler_2018/processed_inputs/
#cp P53_HUMAN_full_11-26-2021_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/P53_HUMAN_Kotler_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/P84126_THETH_Chan_2017/processed_inputs/
#cp P84126_THETH_full_11-26-2021_b04.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/P84126_THETH_Chan_2017/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/PA_I34A1_Wu_2015/processed_inputs/
#cp PA_I34A1_full_theta0.99_04-29-2022_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/PA_I34A1_Wu_2015/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/PABP_YEAST_Melamed_2013/processed_inputs/
#cp PABP_YEAST_full_11-26-2021_b07.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/PABP_YEAST_Melamed_2013/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/POLG_CXB3N_Mattenberger_2021/processed_inputs/
#cp POLG_CXB3N_1-861_theta0.99_04-29-2022_b07.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/POLG_CXB3N_Mattenberger_2021/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/POLG_HCVJF_Qi_2014/processed_inputs/
#cp POLG_HCVJF_theta0.99_1984-2089_11-26-2021_b08.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/POLG_HCVJF_Qi_2014/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/PTEN_HUMAN_Matreyek_2021/processed_inputs/
#cp PTEN_HUMAN_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/PTEN_HUMAN_Matreyek_2021/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/PTEN_HUMAN_Mighell_2018/processed_inputs/
#cp PTEN_HUMAN_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/PTEN_HUMAN_Mighell_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/Q2N0S5_9HIV1_Haddox_2018/processed_inputs/
#cp Q2N0S5_9HIV1_full_theta0.99_04-29-2022_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/Q2N0S5_9HIV1_Haddox_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/Q59976_STRSQ_Romero_2015/processed_inputs/
#cp Q59976_STRSQ_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/Q59976_STRSQ_Romero_2015/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/R1AB_SARS2_Flynn_growth_2022/processed_inputs/
#cp R1AB_SARS2_02-19-2022_b07.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/R1AB_SARS2_Flynn_growth_2022/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/RASH_HUMAN_Bandaru_2017/processed_inputs/
#cp RASH_HUMAN_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/RASH_HUMAN_Bandaru_2017/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/REV_HV1H2_Fernandes_2016/processed_inputs/
#cp REV_HV1H2_full_theta0.99_04-29-2022_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/REV_HV1H2_Fernandes_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/RL401_YEAST_Mavor_2016/processed_inputs/
#cp RL401_YEAST_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/RL401_YEAST_Mavor_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/RL401_YEAST_Roscoe_2013/processed_inputs/
#cp RL401_YEAST_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/RL401_YEAST_Roscoe_2013/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/RL401_YEAST_Roscoe_2014/processed_inputs/
#cp RL401_YEAST_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/RL401_YEAST_Roscoe_2014/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SC6A4_HUMAN_Young_2021/processed_inputs/
#cp SC6A4_HUMAN_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SC6A4_HUMAN_Young_2021/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SCN5A_HUMAN_Glazer_2019/processed_inputs/
#cp SCN5A_HUMAN_1611-1642_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SCN5A_HUMAN_Glazer_2019/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SPG1_STRSG_Olson_2014/processed_inputs/
#cp SPG1_STRSG_full_11-26-2021_b07.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SPG1_STRSG_Olson_2014/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SPIKE_SARS2_Starr_bind_2020/processed_inputs/
#cp SPIKE_SARS2_theta0.99_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SPIKE_SARS2_Starr_bind_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SPIKE_SARS2_Starr_expr_2020/processed_inputs/
#cp SPIKE_SARS2_theta0.99_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SPIKE_SARS2_Starr_expr_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SRC_HUMAN_Ahler_CD_2019/processed_inputs/
#cp SRC_HUMAN_full_11-26-2021_b06.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SRC_HUMAN_Ahler_CD_2019/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SUMO1_HUMAN_Weile_2017/processed_inputs/
#cp SUMO1_HUMAN_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SUMO1_HUMAN_Weile_2017/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SYUA_HUMAN_Newberry_2020/processed_inputs/
#cp SYUA_HUMAN_full_04-29-2022_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/SYUA_HUMAN_Newberry_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TADBP_HUMAN_Bolognesi_2019/processed_inputs/
#cp TADBP_HUMAN_full_11-26-2021_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TADBP_HUMAN_Bolognesi_2019/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TAT_HV1BR_Fernandes_2016/processed_inputs/
#cp TAT_HV1BR_full_theta0.99_04-29-2022_b09.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TAT_HV1BR_Fernandes_2016/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TPK1_HUMAN_Weile_2017/processed_inputs/
#cp TPK1_HUMAN_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TPK1_HUMAN_Weile_2017/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TPMT_HUMAN_Matreyek_2018/processed_inputs/
#cp TPMT_HUMAN_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TPMT_HUMAN_Matreyek_2018/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TPOR_HUMAN_Bridgford_S505N_2020/processed_inputs/
#cp TPOR_HUMAN_full_11-26-2021_b01.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TPOR_HUMAN_Bridgford_S505N_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TRPC_SACS2_Chan_2017/processed_inputs/
#cp TRPC_SACS2_full_11-26-2021_b07.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TRPC_SACS2_Chan_2017/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TRPC_THEMA_Chan_2017/processed_inputs/
#cp TRPC_THEMA_full_11-26-2021_b07.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/TRPC_THEMA_Chan_2017/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/UBC9_HUMAN_Weile_2017/processed_inputs/
#cp UBC9_HUMAN_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/UBC9_HUMAN_Weile_2017/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/UBE4B_MOUSE_Starita_2013/processed_inputs/
#cp UBE4B_MOUSE_full_11-26-2021_b05.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/UBE4B_MOUSE_Starita_2013/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/VKOR1_HUMAN_Chiasson_abundance_2020/processed_inputs/
#cp VKOR1_HUMAN_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/VKOR1_HUMAN_Chiasson_abundance_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/VKOR1_HUMAN_Chiasson_activity_2020/processed_inputs/
#cp VKOR1_HUMAN_full_11-26-2021_b03.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/VKOR1_HUMAN_Chiasson_activity_2020/processed_inputs/nt_alignment_msa.fna
#
#mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/YAP1_HUMAN_Araya_2012/processed_inputs/
#cp YAP1_HUMAN_full_11-26-2021_b02.a2m /groups/doudna/projects/daniel_projects/prywes_n/proteingym_dms_0.70/YAP1_HUMAN_Araya_2012/processed_inputs/nt_alignment_msa.fna

