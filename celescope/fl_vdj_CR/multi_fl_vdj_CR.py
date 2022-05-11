from celescope.tools.multi import Multi
from celescope.fl_vdj_CR.__init__ import __ASSAY__

class Multi_fl_vdj_CR(Multi):
    """

    ## Download and unpack cellranger soft and reference file.
    ```
    wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-vdj/cellranger-6.1.2.tar.gz?Expires=1646072261&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC12ZGovY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NDYwNzIyNjF9fX1dfQ__&Signature=Z-2m906CV5Rb1snIAga-QDSXYSZ8cNqCj1EECGP4uloU3qH~uCMH42MHf4TNnDL2zAsKA7cXsCsQYz0A9yJdNh7dfRT8ohpuAzASFx5Pj-bkqfw4p2tql55IIaPN0zqxyUuyZ9sfKl5qTQX82LoVolRpiBUL8dF9nr~bA2P1gJZ~xg1QssS7icR5MmTzvKKS5NYkezG8vWaTiEdXU0nuKI2ciZSX5GOMeIRW-YYR7mJwHmBbTVxe0o-uBuUtqor0Y98jdIv8Z~dwMjujRjrEShdCGNixTSonGzeS2~9CXqWquCJIOolqFFkcFHgXkD7ZWNfSXWbTxuF57rCsub98pA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

    Reference: human and mouse
    wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
    wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0.tar.gz

    tar -xzvf cellranger-6.1.2.tar.gz
    tar -xzvf refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
    tar -xzvf refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0.tar.gz
    ```

    ## Usage
    
    ```
    conda activate celescope
    multi_fl_vdj_CR \\
        --mapfile ./test.mapfile \\
        --chemistry flv \\
        --mem 10 \\
        --thread 8 \\
        --allowNoLinker \\
        --seqtype TCR \\
        --ref_path "/SGRNJ/Database/script/soft/cellranger/vdj_ref/6.0.0/hs/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0" \\
        --soft_path "/SGRNJ03/randd/cjj/soft/cellranger/cellranger-6.1.2/cellranger" 
    ```
    """

    def convert(self, sample):
        step = 'convert'
        cmd_line = self.get_cmd_line(step, sample)
        out_dir = f'{self.outdir_dic[sample][step]}'
        fq2 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq'
        cmd = (
            f'{cmd_line} '
            f'--fq2 {fq2} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=1)
        self.process_cmd(cmd, step, sample, m=5, x=1)


    def assemble(self, sample):
        step = 'assemble'
        cmd_line = self.get_cmd_line(step, sample)
        out_dir = f'{self.outdir_dic[sample][step]}'
        fqs_dir = f'{self.outdir_dic[sample]["convert"]}'
        cmd = (
            f'{cmd_line} '
            f'--fqs_dir {fqs_dir} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=self.args.thread)
        self.process_cmd(cmd, step, sample, m=self.args.mem, x=self.args.thread)


    def annotation(self,sample):
        step = 'annotation'
        cmd_line = self.get_cmd_line(step, sample)
        out_dir = f'{self.outdir_dic[sample][step]}'
        barcode_dict = f'{self.outdir_dic[sample]["convert"]}/barcode_correspond.txt'
        cmd=(
            f'{cmd_line} '
            f'--barcode_dict {barcode_dict} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=self.args.thread)
        self.process_cmd(cmd, step, sample, m=8, x=self.args.thread)


    def match(self,sample):
        step = 'match'
        cmd_line = self.get_cmd_line(step, sample)
        out_dir = f'{self.outdir_dic[sample][step]}'
        barcode_dict = f'{self.outdir_dic[sample]["convert"]}/barcode_correspond.txt'
        match_dir = f'{self.col4_dict[sample]}'
        cmd = (
            f'{cmd_line} '
            f'--barcode_dict {barcode_dict} '
            f'--match_dir {match_dir} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=self.args.thread)
        self.process_cmd(cmd, step, sample, m=8, x= self.args.thread)
    

    def summarize(self, sample):
        step = 'summarize'
        cmd_line = self.get_cmd_line(step, sample)
        out_dir = f'{self.outdir_dic[sample][step]}'
        barcode_dict = f'{self.outdir_dic[sample]["convert"]}/barcode_correspond.txt'
        cmd=(
            f'{cmd_line} '
            f'--barcode_dict {barcode_dict} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=self.args.thread)
        self.process_cmd(cmd, step, sample, m=8, x=self.args.thread)
    
    
    def mapping(self,sample):
        step = 'mapping'
        cmd_line = self.get_cmd_line(step,sample)
        out_dir = f'{self.outdir_dic[sample][step]}'
        match_dir = f'{self.col4_dict[sample]}'
        cmd=(
            f'{cmd_line} '
            f'--match_dir {match_dir} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=1)
        self.process_cmd(cmd, step, sample, m=5, x=1)

def main():
    multi = Multi_fl_vdj_CR(__ASSAY__)
    multi.run()
    

if __name__ == '__main__':
    main()
