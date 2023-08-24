## Filogenetyczny pipeline
Konstrukcja drzew gatunków i superdrzew dla zadanej listy numerów dostępu proteomów z NCBI. Białka są klastrowane, klastry poddawane uliniowieniu (przy zastosowaniu parametru --bijective również modyfikacji, która poprzez porównania sekwencji w obrębie klastra usuwa pewne sekwencje, jeśli występują duplikaty z jednego gatunku), a następnie metodą Nejghbour Joining dla każdego klastra konsturowane jest drzewo genów. Z usyskanych drzew obliczane są drzewa konsensusowe i super drzewa.

### Dependencies

* MMseqs2 (https://github.com/soedinglab/MMseqs2)
* MAFFT (https://mafft.cbrc.jp/alignment/software/linux.html)
* Deduptree (https://genome.cs.iastate.edu/DupTree)
* NCBI datasets (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)

### Skrypty

* #### download_prot_acc.sh 
pobiera proteomy, których numery dostępu wypisane są w kolejnych linijkach pliku accession_nums.txt - plik ten musi się znajdować w tym samym katalogu co skypt

* #### clusters_to_gene_trees.py
Przyjmuje output mmsesq2 (plik fasta zawierający wszystkie sekwencje) i zwraca uliniowienia klastrów i drzewa genów. Jako parametry przyjmuje:
  1. Ścieżkę do pliku z klastrami (output mmseqs2)
  2. Ścieżkę do katalogu, w którym zapisywane mają być wyniki
  3. Minimalna liczba genomów, które mają występować w klastrze, aby był analizowany.
  4. Parametr --bijective, który oznacza, że klastry mają zostać poddane obróbce tak, aby w każdym występowała maksymalnie jedna sekwencja z każdego genomu.

* #### trees_to_one_file.py
Jako parametr przyjmuje ścieżkę do katalogu, w którym znajdują się tylko drzewa genów oraz ścieżkę do katalogu do którego mają zostać zapisane wyniki. Zapisuje wszystkie drzewa do jednego pliku, dodatkowo odfiltrowując te, z ujemnymi długościami krawędzi i poddając nazwy genomów obróbce, aby plik był akceptowany przez duptree.

* #### consensus_trees.py
Jako parametry przyjmuje ścieżkę do pliku newick z drzewami i ścieżkę do katalogu, do którego mają zostać zapisane wyniki (drzewa konsensusowe).

* #### constr_consensus.r
Konstruuje drzewo konsensusowe dla pliku newick z drzewami.

* #### run_pipeline.sh
Uruchamia cały pipeline, wynikiem są: drzewo konsensusowe i dwa super drzewa. Musi być w tym samym katalogu co wszystkie wymienione wyżej skrypty, oraz duptree i datasets.
