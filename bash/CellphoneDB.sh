## source python virtualenv
source $HOME/my_python-3.6.3/bin/activate


## do run
cellphonedb method statistical_analysis /~/data/seurat/cellphonedb_meta.txt /~/data/seurat/cellphonedb_count.txt --counts-data=gene_name --project-name=IRI --iterations=1000 --threshold=0.05 --output-path=/~/CellphoneDB/ --threads=14 --result-precision=4