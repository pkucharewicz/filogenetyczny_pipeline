mkdir proteoms_download
cd proteoms_download

while IFS="" read -r name || [ -n "$p" ]
do
	/home/pampuch/studia/genomika_porównawcza/GP_projekt/final_pipeline/datasets download genome accession "$name" --filename "$name.zip" --exclude-rna --exclude-genomic-cds --exclude-gff3 --exclude-seq
	
	mkdir temp_prot
	mv "$name.zip" temp_prot
	cd temp_prot
	unzip "$name.zip" 
	mv ncbi_dataset/data/GCA*/protein.faa ../${name// /_}.faa
	cd ..
	rm -r temp_prot
	
done < ../accession_nums.txt

cat GCA* > all_proteoms.fasta
mv all_proteoms.fasta ..


#name="Thermosynechococcus elongatus BP-1"
#/home/pampuch/studia/genomika_porównawcza/GP_projekt/datasets download genome taxon "$name" --filename "${name// /_}.zip" --exclude-rna --exclude-genomic-cds --exclude-gff3 --exclude-seq --assembly-source 'GenBank'

#unzip "$name.zip"

#mv ncbi_dataset/data/GCA*/protein.faa .
#rm -r ncbi_dataset
#mv protein.faa ${name// /_}.faa