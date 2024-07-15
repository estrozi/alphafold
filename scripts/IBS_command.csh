#!/bin/tcsh -f
# Author: Leandro F. Estrozi, Institut de Biologie Structurale, Grenoble, CNRS.
# This script runs Alphafold2 and speed-up things by taking into account that:
# 1) Many AF2 predictions are already available in public databases, so it queries them to check that first.
# 2) MSA calculations can be done only once for any given sequence, thus a copy of historical/previous results
# are stored in /storage/Data/AF2seqs.txt and /storage/Data/AF2msas/tmp.??????????
# 3) Because there are 5 AI-models, they can be run in parallel if they fit in the VRAM.
#
# This script places its outputs in the current folder but it also uses /storage/Data/ as a temporary space.

  if($#argv < 2) then
echo "Usage: ${0:t} fastafile #predictions_per_model [seed] [--skipqueryDBs] [--ignoretags] [--notemplates]";
echo "";
echo "In the AF2 context, a 'model' is NOT a PDB/structure. A model is a neural";
echo "network with many layers, which is able to learn from data to do specific";
echo "tasks. AF2 has 5 slightly different *models* to make structure predictions.";
echo "";
echo "So here, a PDB/strucutre will be called a PREDICTION, not a 'model'";
echo "";
echo "The choice of the number of predictions per model depends on the stability";
echo "of the obtained results. When the MSAs are shalow (low homology) or if the";
echo "predictions come out with low confidence, a bigger (>2) number of predictions";
echo "per model becomes necessary.";
echo "On the contrary, when there are lots of homology and the predictions come";
echo "with high confidence, no more than 2 predictions per model are needed.";
echo "";
echo "Keep in mind that the excution time depends (almost) linearly on this number";
echo "so values like 5 or more often make calculations pretty long for no reason.";
echo "I recommend to begin with 2 predictions per model, inspect the results and";
echo "only later go for more predictions if needed.";
exit 0;
  endif

set msa_dbs_mode = full_dbs;
#set msa_dbs_mode = reduced_dbs;

set fasta = $1;
  if( ! -e $fasta || -z $fasta ) then
echo "ERROR: fastafile not found";
exit 1;
  endif

#The sbgrid option below depends on 'mount --bind' commands to make the adapted AF2 scripts visible.
# So root privileges are necessary.
#The docker option tries to avoid that. For the docker option you need first:
# 1) To clone the IBS AF2 repository:
# #>git clone https://github.com/estrozi/alphafold.git
# 2) To build the IBS AF2 Docker image:
# #>cd alphafold/
# #>docker build --build-arg UID=`id -u` --build-arg GID=`id -g` --build-arg UNAME=`id -un` -f docker/Dockerfile -t alphafold .
# 3) To test if it has access to the GPUs:
# #>docker run --rm --gpus all nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu18.04 nvidia-smi
#
# If you get a normal nvidia-smi output, you probably have eveything you need.
#set INSTALLATION_TYPE = "sbgrid";
set INSTALLATION_TYPE = "docker";

set AF2_CMD = "run_alphafold.py";
set prefix = "";

  if( $INSTALLATION_TYPE == "docker" ) then
# not sure the -v blabla socket: is necessary
#set AF2_CMD = "docker run --user "`id -un`":"`id -gn`" --rm --gpus all --mount type=bind,source=/storage,target=/app/storage -v /var/run/nslcd/socket:/var/run/nslcd/socket";
set AF2_CMD = "docker run --user "`id -un`":"`id -gn`" --rm --gpus all --mount type=bind,source=/storage,target=/app/storage";
set prefix = "/app";
  endif

set n_preds = $2;

set invalid_letters = `grep -v '^>' $fasta | tr '[:lower:]' '[:upper:]' | grep -e 'J' -e 'O' |  wc -l`;
if($status) exit 1;
  if( $invalid_letters ) then
#maybe other letters are also considered invalid by AF2
echo "ERROR: there are invalid characters (J,O) in the sequence(s)";
exit 1;
  endif

set len = `grep -v '^>' $fasta | tr -d '\n' | sed -e 's/\s//g' | wc -m`;
if($status) exit 1;
set nsq = `grep -v '^>' $fasta | wc -l`;
if($status) exit 1;

awk 'function ltrim(s){sub(/^[ \t\r\n]+/,"",s);return s}function rtrim(s){sub(/[ \t\r\n]+$/,"",s);return s}function trim(s){return rtrim(ltrim(s))}BEGIN{seq=""}{if(substr($1,1,1)==">"){if(seq!=""){aarray[seq]++};seq=""} else {seq=seq""toupper(trim($0))}}END{aarray[seq]++;for(seq in aarray){print seq}}' $fasta >! $$.unique_seqs;
if($status) exit 1;
set nus = `cat $$.unique_seqs | wc -l`;
if($status) exit 1;

set skipqueryDBs = 0;
  if($#argv > 3) then
    if($4 == "--skipqueryDBs") then
set skipqueryDBs = 1;
echo "Skipping query DBs";
    endif
  endif

set ignoretags = 0;
  if($#argv > 4) then
    if("$5" == "--ignoretags") then
set ignoretags = 1;
    endif
  endif

  if( ! $ignoretags ) then
#Search for tags 
set i = 1;
set tagname = ( his-tag HQ-tag    HN-tag        HAT-tag           SpyTag        KTag      SnoopTag     SpyTag002     SnoopTagJr          DogTag             AVI-atg );
    foreach tag ( HHHHHH HQHQHQ HNHNHNHNHNHN KDHLIHNVHKEEHAHAHNK AHIVMVDAYKPTK ATHIKFSKRD KLGDIEFIKVNK VPTIVMVDAYKRYK KLGSIEFIKVNK DIPATYEFTDGKHYITNEPIPPK GLNDIFEAQKIEWHE )
grep -n $tag $$.unique_seqs;
      if( $status == 0 ) then
echo "The sequence(s) above seem(s) to have a $tagname[$i] $tag";
echo "Please verify that is what you want and run with the flag --ignoretags";
exit 1;
      endif
@ i++;
    end
  endif

set outdir = ${fasta:t:r}_$msa_dbs_mode;
set outdir2 = /storage/Data/$outdir;
  if( -e $outdir) then
    if( ! -e $outdir2 ) then
echo "Moving already existing output folder $outdir to $outdir2"
echo "(it will be moved back in the end if everything goes fine)"
mv $outdir $outdir2;
if($status) exit 1;
    else
echo "Conflict: Both $outdir and $outdir2 exit. Aborting...";
exit 1;
    endif
  else
mkdir $outdir2;
  endif
if($status) exit 1;

  if( -e $outdir2/${fasta:t:r}/ranking_debug.json && ! -z $outdir2/${fasta:t:r}/ranking_debug.json ) then
echo "Final file ranking_debug.json found. Aborting..."
echo "Moving $outdir2 to $outdir"
mv $outdir2 $outdir;
if($status) exit 1;
exit 1;
  endif

set alphabet = ( A B C D E F G H I J K L M N O P Q R S T U V W X Y Z );

#AF2seqs.txt has 1st column: sequence 2nd column: msa_folder/ with the .sto/.a3m files inside.
  if (! -e /storage/Data/AF2seqs.txt) then
touch /storage/Data/AF2seqs.txt;
if($status) exit 1;
chmod 664 /storage/Data/AF2seqs.txt;
if($status) exit 1;
  else
set i = 1;
if($status) exit 1;
    foreach seq ( `cat $$.unique_seqs` )
      if( `awk -v seq=$seq 'BEGIN{p=0}{if(seq==$1) {print 1; p=1; exit;}}END{if(p==0) {print 0} }' /storage/Data/AF2seqs.txt` ) then
mkdir -p $outdir2/${fasta:t:r}/msas/$alphabet[$i]/;
cp -i -a `awk -v seq=$seq '{if(seq==$1) {print $2"/*"}}' /storage/Data/AF2seqs.txt` $outdir2/${fasta:t:r}/msas/$alphabet[$i]/
if($status) exit 1;
      endif
@ i++
    end
  endif

#all possible templates
set templatedate = 3000-01-01
  if($#argv > 5) then
    if("$6" == "--notemplates") then
#no templates
set templatedate = 0001-01-01
    else
echo "${0:t}: ERROR: 5th argument is not --notemplates";
exit 1;
    endif
  endif

set seed = "";
  if($#argv > 2) then
    if( $3 != "") then
      if(`echo $3 | awk '{if(math($1,/[^0-9]+/)){print 1} else {print 0}}'`) then
echo "${0:t}: ERROR 3rd argument is not a numeric seed";
exit 1;
      else
set seed = "--random_seed $3";
      endif
    endif
  endif

#relax is not working with the current setup (Debian 12 + nvidia drivers + cuda version)
#disabled for now
set common_args = "$seed --models_to_relax none --num_multimer_predictions_per_model $n_preds"
set data = $prefix/storage/Alphafold/data;

# preset monomer case
  if (0) then
set common_args = "$seed --models_to_relax none --model_preset=monomer"
set preset_dependent_args = "--pdb70_database_path $data/pdb70/pdb70"
  else
set preset_dependent_args = "--pdb_seqres_database_path $data/pdb_seqres/pdb_seqres.txt --uniprot_database_path $data/uniprot/uniprot.fasta"
  endif
  if ( $msa_dbs_mode == "full_dbs" ) then
set msa_dbs_mode_dependent_args = "$common_args --bfd_database_path $data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt --uniref30_database_path $data/uniref30/UniRef30_2021_03"
  else
set msa_dbs_mode_dependent_args = "$common_args --small_bfd_database_path $data/small_bfd/bfd-first_non_consensus_sequences.fasta"
  endif
#agressive memory use
#setenv XLA_PYTHON_CLIENT_PREALLOCATE true
##setenv XLA_PYTHON_CLIENT_MEM_FRACTION 75.0 #default
#setenv XLA_PYTHON_CLIENT_MEM_FRACTION 100.0
#conservative memory use
#setenv XLA_PYTHON_CLIENT_PREALLOCATE false
#setenv TF_FORCE_UNIFIED_MEMORY 1
#setenv XLA_PYTHON_CLIENT_MEM_FRACTION 0.75
#setenv XLA_PYTHON_CLIENT_ALLOCATOR platform
##setenv CUDA_VISIBLE_DEVICES 0

if ($skipqueryDBs) goto labelskipqueryDBs;

set i = 1;
set ids = "";
set nhits = 0;
  while ($i <= $nus)
set cur_seq = `head -$i $$.unique_seqs | tail -1`;
#the commands below should be able to verify if a sequence is already present in the PDB
echo "Querying PDB server for EXPERIMENTAL sequence ${cur_seq}";
curl -s X GET --header 'Accept:application/json' 'https://search.rcsb.org/rcsbsearch/v2/query?json=%7B%0A%20%20%22query%22%3A%20%7B%0A%20%20%20%20%22type%22%3A%20%22group%22%2C%0A%20%20%20%20%22logical_operator%22%3A%20%22and%22%2C%0A%20%20%20%20%22nodes%22%3A%20%5B%20%7B%0A%20%20%20%20%20%20%22type%22%3A%22terminal%22%2C%0A%20%20%20%20%20%20%22service%22%3A%22seqmotif%22%2C%0A%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%22value%22%3A%22'${cur_seq}'%22%2C%0A%20%20%20%20%20%20%20%20%22pattern_type%22%3A%22simple%22%2C%0A%20%20%20%20%20%20%20%20%22sequence_type%22%3A%22protein%22%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%7D%20%5D%0A%20%20%7D%2C%0A%20%20%22return_type%22%3A%20%22polymer_entity%22%2C%0A%20%20%22request_options%22%3A%20%7B%0A%20%20%20%20%22results_content_type%22%3A%20%5B%22experimental%22%5D%2C%0A%20%20%20%20%22scoring_strategy%22%3A%20%22combined%22%0A%20%20%7D%0A%7D' | awk -v FS=',' '{for(i=1; i<=NF; i++) print $i}' >! $$.fetchPDB.$i;
    if( ! -z $$.fetchPDB.$i ) then
grep 'identifier' $$.fetchPDB.$i | sed -e 's/"result_set":\[//' >! $$.fetchPDB_ids.$i;
grep 'score' $$.fetchPDB.$i >! $$.fetchPDB_scores.$i;
      if( ! -z $$.fetchPDB_ids.$i && ! -z $$.fetchPDB_scores.$i ) then
paste -d ' ' $$.fetchPDB_ids.$i $$.fetchPDB_scores.$i | sed -e 's/:/ /g' >! $$.fetchPDB.$i;
        if( `sort -n -k 4 $$.fetchPDB.$i | tail -1 | awk '{if($4>=0.999){print 1} else {print 0}}'` ) then
@ nhits++
echo "Sequence $i match with PDB https://www.rcsb.org/structure/"`sort -n -k 4 $$.fetchPDB.$i | tail -1 | awk '{print substr($2,2,4)}'`;
echo -n " $i" >>! $$.matchPDB.`sort -n -k 4 $$.fetchPDB.$i | tail -1 | awk '{print substr($2,2,4)}'`;
set ids = "$ids "$i;
        endif
      endif
rm -f $$.fetchPDB.$i $$.fetchPDB_ids.$i $$.fetchPDB_scores.$i;
    else
rm -f $$.fetchPDB.$i;
    endif
@ i++
  end
  if( "$ids" != "" ) then
    foreach f ( $$.matchPDB.???? )
echo "" >> $f;
      if ( "`cat $f`" == "$ids" ) then
echo "All the requested sequences match the same PDB: https://www.rcsb.org/structure/${f:e} Aborting...";
rm -f $$.matchPDB.????;
exit 1;
      endif
    end
  endif

rm -f $$.matchPDB.????;

set i = 1;
set ids = "";
set nhits = 0;
  while ($i <= $nus)
set cur_seq = `head -$i $$.unique_seqs | tail -1`;
#the commands below should be able to verify if a sequence is already present in the AFDB
echo "Querying PDB server for COMPUTATIONAL sequence ${cur_seq}";
curl -s X GET --header 'Accept:application/json' 'https://search.rcsb.org/rcsbsearch/v2/query?json=%7B%0A%20%20%22query%22%3A%20%7B%0A%20%20%20%20%22type%22%3A%20%22group%22%2C%0A%20%20%20%20%22logical_operator%22%3A%20%22and%22%2C%0A%20%20%20%20%22nodes%22%3A%20%5B%20%7B%0A%20%20%20%20%20%20%22type%22%3A%22terminal%22%2C%0A%20%20%20%20%20%20%22service%22%3A%22seqmotif%22%2C%0A%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%22value%22%3A%22'${cur_seq}'%22%2C%0A%20%20%20%20%20%20%20%20%22pattern_type%22%3A%22simple%22%2C%0A%20%20%20%20%20%20%20%20%22sequence_type%22%3A%22protein%22%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%7D%20%5D%0A%20%20%7D%2C%0A%20%20%22return_type%22%3A%20%22polymer_entity%22%2C%0A%20%20%22request_options%22%3A%20%7B%0A%20%20%20%20%22results_content_type%22%3A%20%5B%22computational%22%5D%2C%0A%20%20%20%20%22scoring_strategy%22%3A%20%22combined%22%0A%20%20%7D%0A%7D' | awk -v FS=',' '{for(i=1; i<=NF; i++) print $i}' >! $$.fetchAFDB.$i;
    if( ! -z $$.fetchAFDB.$i ) then
grep 'identifier' $$.fetchAFDB.$i | sed -e 's/"result_set":\[//' >! $$.fetchAFDB_ids.$i;
grep 'score' $$.fetchAFDB.$i >! $$.fetchAFDB_scores.$i;
      if( ! -z $$.fetchAFDB_ids.$i && ! -z $$.fetchAFDB_scores.$i ) then
paste -d ' ' $$.fetchAFDB_ids.$i $$.fetchAFDB_scores.$i | sed -e 's/:/ /g' >! $$.fetchAFDB.$i;
        if( `sort -n -k 4 $$.fetchAFDB.$i | tail -1 | awk '{if($4>=0.999){print 1} else {print 0}}'` ) then
@ nhits++
echo "Sequence $i match with AFDB https://alphafold.ebi.ac.uk/entry/"`sort -n -k 4 $$.fetchAFDB.$i | tail -1 | awk '{print substr($2,7,10)}'`;
echo -n " $i" >>! $$.matchAFDB.`sort -n -k 4 $$.fetchAFDB.$i | tail -1 | awk '{print substr($2,7,10)}'`;
set ids = "$ids "$i;
        endif
      endif
rm -f $$.fetchAFDB.$i $$.fetchAFDB_ids.$i $$.fetchAFDB_scores.$i;
    else
rm -f $$.fetchAFDB.$i;
    endif
@ i++
  end
  if( "$ids" != "" ) then
    foreach f ( $$.matchAFDB.?????????? )
echo "" >> $f;
      if ( "`cat $f`" == "$ids" ) then
echo "All the requested sequences match the same AFDB: https://alphafold.ebi.ac.uk/entry/${f:e} Aborting...";
rm -f $$.matchAFDB.??????????;
exit 1;
      endif
    end
  endif

rm -f $$.matchAFDB.??????????;

labelskipqueryDBs:

cp -i -a $fasta /storage/Data/${fasta:t};
# the -i here is because perhaps the presence of the /storage/Data/blabla.fasta can be the reason to skip

set i = 1;
set nhits = 0;
echo "Job length: $len; #Sequences: $nsq; #UniqueSequences: $nus;"
#optimized for parallel models by LFE (stopat 0 use no GPU):
setenv XLA_PYTHON_CLIENT_PREALLOCATE false
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 1.0
  if( $INSTALLATION_TYPE == "docker" ) then
eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=1.0 alphafold --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat 0" |& tee -a $outdir.AF2_IBS.log;
if($status) exit 1;
  else
$AF2_CMD --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat 0 |& tee -a $outdir.AF2_IBS.log;
if($status) exit 1;
  endif
    
  if($len < 500) then
#try 5 processes in parallel
#optimized for parallel models by LFE:
setenv XLA_PYTHON_CLIENT_PREALLOCATE false
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 20.0
set i = 1;
    while ( $i < 6 )
      if( $INSTALLATION_TYPE == "docker" ) then
( eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=20.0 alphafold --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" & ) |& tee -a $outdir.AF2_IBS.$i.log &;
      else
( $AF2_CMD --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i & ) |& tee -a $outdir.AF2_IBS.$i.log &;
      endif
@ i++
    end
wait;
  else
    if($len < 600) then
#try 3+2 processes in parallel
setenv XLA_PYTHON_CLIENT_PREALLOCATE false
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 33.333
set i = 1;
      while ( $i < 4 )
        if( $INSTALLATION_TYPE == "docker" ) then
( eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=33.333 alphafold --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" & ) |& tee -a $outdir.AF2_IBS.$i.log &;
       else
( $AF2_CMD --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i & ) |& tee -a $outdir.AF2_IBS.$i.log &;
       endif
@ i++
      end
wait;
#try 2 processes in parallel
setenv XLA_PYTHON_CLIENT_PREALLOCATE false
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 50.0
      while ( $i < 6 )
        if( $INSTALLATION_TYPE == "docker" ) then
( eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=50.0 alphafold --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" & ) |& tee -a $outdir.AF2_IBS.$i.log &;
        else
( $AF2_CMD --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i & ) |& tee -a $outdir.AF2_IBS.$i.log &;
        endif
@ i++
      end
wait;
    else
      if($len < 700) then
#try 2+2+1 processes in parallel
#optimized for parallel models by LFE:
setenv XLA_PYTHON_CLIENT_PREALLOCATE false
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 50.0
set i = 1;
        while ( $i < 3 )
          if( $INSTALLATION_TYPE == "docker" ) then
( eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=50.0 alphafold --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" & ) |& tee -a $outdir.AF2_IBS.$i.log &;
          else
( $AF2_CMD --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i & ) |& tee -a $outdir.AF2_IBS.$i.log &;
          endif
@ i++
        end
wait;
        while ( $i < 5 )
          if( $INSTALLATION_TYPE == "docker" ) then
( eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=50.0 alphafold --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" & ) |& tee -a $outdir.AF2_IBS.$i.log &;
          else
( $AF2_CMD --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i & ) |& tee -a $outdir.AF2_IBS.$i.log &;
          endif
@ i++
        end
wait;
#agressive memory use
setenv XLA_PYTHON_CLIENT_PREALLOCATE true
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 100.0
        if( $INSTALLATION_TYPE == "docker" ) then
eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=true --env XLA_PYTHON_CLIENT_MEM_FRACTION=100.0 alphafold --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" |& tee -a $outdir.AF2_IBS.$i.log;
        else
$AF2_CMD --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i |& tee -a $outdir.AF2_IBS.$i.log;
        endif
      endif
    endif
  endif
  if($len < 700) then
set i = 1;
    while ( $i < 6 )
cat $outdir.AF2_IBS.$i.log >>! $outdir.AF2_IBS.log;
rm -f $outdir.AF2_IBS.$i.log;
@ i++;
    end
  endif
#unexpected good thing: if any of the run-in-background inferences above fail because lack of vRAM
#the run_alphafold below run the missing ones one-by-one
#agressive memory use
setenv XLA_PYTHON_CLIENT_PREALLOCATE true
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 100.0
  if( $INSTALLATION_TYPE == "docker" ) then
eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=true --env XLA_PYTHON_CLIENT_MEM_FRACTION=100.0 alphafold --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat 6" |& tee -a $outdir.AF2_IBS.log;
  else
$AF2_CMD --fasta_paths $prefix/storage/Data/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat 6 |& tee -a $outdir.AF2_IBS.log;
  endif
if($status) exit 1;

set i = 1;
if($status) exit 1;
  foreach seq ( `cat $$.unique_seqs` )
grep -i "$seq" /storage/Data/AF2seqs.txt >& /dev/null;
    if($status == 1) then #match not found
set tmpdir = `mktemp -u -d -t -p /storage/Data/AF2msas/`;
if($status) exit 1;
mkdir -p $tmpdir;
if($status) exit 1;
chmod 775 $tmpdir;
if($status) exit 1;
echo "$seq $tmpdir" >>! /storage/Data/AF2seqs.txt;
if($status) exit 1;
cp -a $outdir2/${fasta:t:r}/msas/$alphabet[$i]/* $tmpdir/; 
if($status) exit 1;
    endif
@ i++
  end

echo "Moving $outdir2 to $outdir"
mv $outdir2 $outdir;
if($status) exit 1;
rm -f /storage/Data/${fasta:t};
if($status) exit 1;
rm -f $$.unique_seqs;
if($status) exit 1;

exit 0;
