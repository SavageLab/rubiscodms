# == Native Modules ==
# == Installed Modules ==
import pandas as pd
# == Project Modules ==


def main():
    # === SNAKEMAKE I/O ===
    # == INPUTS
    barcode_report_list = list(snakemake.input.barcode_report_list)
    # == OUTPUTS
    concat_counts_path = str(snakemake.output.concat_counts_path)
    # DEBUG
    # barcode_report_list = ["/groups/doudna/projects/daniel_projects/rubiscodms/rubisco_nextseq_processing/barcodes/NP_11_66_8_pacbioMerged_barcodeCounts.csv",
    #                        "/groups/doudna/projects/daniel_projects/rubiscodms/rubisco_nextseq_processing/barcodes/NP_11_66_15_pacbioMerged_barcodeCounts.csv"]

    loop_report_entry = pd.read_csv(barcode_report_list[0],index_col=0)
    for list_index in range(1, len(barcode_report_list)):
        df_pacbio_barcode = pd.read_csv(barcode_report_list[list_index], index_col=0)
        loop_report_entry = pd.concat([df_pacbio_barcode, loop_report_entry])
    df_concat_barcode = loop_report_entry.convert_dtypes()

    df_concat_barcode.to_csv(concat_counts_path)


if __name__ == "__main__":
    main()
