# == Native Modules ==
# == Installed Modules ==
import panda as pd
# == Project Modules ==


def main():
    # === SNAKEMAKE I/O ===
    # == INPUTS
    barcode_report_list = str(snakemake.input.barcode_report_list)
    # == PARAMS
    flanking_sequence = str(snakemake.params.flanking_sequence)
    # == OUTPUTS
    concat_counts_path = str(snakemake.output.pacbio_merged_counts)


    loop_report_entry = barcode_report_list[0]
    for list_index in range(1, len(barcode_report_list)):
        df_pacbio_barcode = pd.read_csv(barcode_report_list[list_index], index_col=0)
        loop_report_entry = pd.concat([df_pacbio_barcode, loop_report_entry])
    df_concat_barcode = loop_report_entry.convert_dtypes()



if __name__ == "__main__":
	main()
