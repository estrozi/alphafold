#!/bin/tcsh -f
# Author: Leandro F. Estrozi, Institut de Biologie Structurale, Grenoble, CNRS.
# This script runs Alphafold2 and speed-up things by taking into account that:
# 1) Many AF2 predictions are already available in public databases, so it queries them to check that first.
# 2) MSA calculations can be done only once for any given sequence, thus a copy of historical/previous results
# are stored in /storage/Data/AF2jobs/AF2seqs.txt and /storage/Data/AF2jobs/AF2msas/tmp.??????????
# 3) Because there are 5 AI-models, they can be run in parallel if they fit in the VRAM.
#
# This script places its outputs in the current folder but it also uses /storage/Data/AF2jobs/ as a temporary space.

#REMINDER: the caller add two arguments (uid and jobname) before the fastafile

umask 113; # rw-rw-r--

set caller = "/storage/Alphafold/scripts/alphafold2_caller.bin"
  if($#argv < 4) then
echo "Usage: $caller fastafile #predictions_per_model [seed] [--skipqueryDBs] [--ignoretags] [--notemplates] [model]";
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
exit 1;
  endif

onintr failed;

set msa_dbs_mode = full_dbs;
#set msa_dbs_mode = reduced_dbs;

set u = $1;
set g = 7182;

set jobname = $2;

set fasta = $3;
set outdir = ${fasta:t:r}_$msa_dbs_mode;
set outdir2 = /storage/Data/AF2jobs/$outdir;
  if( ! -e $fasta || -z $fasta ) then
echo "ERROR: fastafile not found" | tee -a $outdir.log;
goto failed;
  endif

  if( -e $outdir) then
    if( ! -e $outdir2 ) then
echo "Moving already existing output folder $outdir to $outdir2" | tee -a $outdir.log;
echo "(it will be moved back in the end if everything goes fine)" | tee -a $outdir.log;
mv $outdir $outdir2;
if($status) goto failed;
    else
echo "Conflict: Both $outdir and $outdir2 exit. Aborting..." | tee -a $outdir.log;
goto failed;
    endif
rm -f $outdir2/finished.txt;
  else
mkdir -m 775 $outdir2;
chown ${u}:${g} $outdir2;
if($status) goto failed;
  endif

chmod 770 $outdir2;
if($status) goto failed;
chgrp 7182 $outdir2;
if($status) goto failed;
echo 0 >! $outdir2/running.txt;
if($status) goto failed;
chown ${u}:${g} $outdir2/running.txt
if($status) goto failed;

grep -v '^>' $fasta | tr '[:lower:]' '[:upper:]' | grep -e 'J' -e 'O' |  wc -l >& /dev/null;
  if(! $status) then
set invalid_letters = `grep -v '^>' $fasta | tr '[:lower:]' '[:upper:]' | grep -e 'J' -e 'O' |  wc -l`;
    if( $invalid_letters ) then
echo "ERROR: there are $invalid_letters invalid characters (J,O) in the sequence(s)" | tee -a $outdir.log;
goto failed;
    endif
  endif

set len = `grep -v '^>' $fasta | tr -d '\n' | sed -e 's/\s//g' | wc -m`;
  if($status) then
echo "ERROR: failed geting length.";
goto failed;
  endif
  if($len <= 0) then
echo "ERROR: length <= 0";
goto failed;
  endif

#FOR ADMINS (to be done once):
#The sbgrid option below depends on 'mount --bind' commands to make the modified AF2 scripts visible
#inside the sbgrid file tree. These mount --bind commands are in the UGA doc file.
# So, root privileges are necessary if you want to run with the sbgrid option.
#The docker option tries to avoid that. For the docker option you need first:
# 1) To clone the IBS AF2 repository:
# #>git clone https://github.com/estrozi/alphafold.git
# 2) To build the IBS AF2 Docker image (as a ADMIN you need to be in docker group):
# #>cd alphafold/
# #>docker build --build-arg GNAME=ibs-gdt-iaaccess --build-arg GID=7182 --build-arg UNAME=`id -un` -f docker/Dockerfile -t alphafold .
# 3) To test if it has access to the GPUs:
# #>docker run --rm --gpus all nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu18.04 nvidia-smi
# If you get a normal nvidia-smi output, you probably have eveything you need.
# For normal users (not belonging to the docker group), we do the following:
# 4) ...

#set INSTALLATION_TYPE = "sbgrid";
set INSTALLATION_TYPE = "docker";

#REMIND: to set RUN_PARALLEL to 1 with sbgrid you need to do the "mount --bind" business.
#REMIND: if RUN_PARALLEL = 0, the argument "model" has no effect.
set RUN_PARALLEL = 1;

  if($RUN_PARALLEL && $INSTALLATION_TYPE == "sbgrid") then
mount | grep af_bin
    if($status) then
echo "ERROR: trying to run in parallel with sbgrid without the: mount --bind" | tee -a $outdir.log;
goto failed;
    endif
  endif

set AF2_CMD = "srun -J ${jobname} --gres=gpu run_alphafold.py";
set prefix = "";

  if( $INSTALLATION_TYPE == "docker" ) then
set gres = "gpu";
#    if($len < 500) then
#set gres = "shard";
#    else
#      if($len < 600) then
#set gres = "shard:3";
#      else
        if($len < 700) then
set gres = "shard";
        endif
#      endif
#    endif
set AF2_CMD = "srun -J ${jobname} --gres=${gres} docker run --user "${u}":"${g}" --rm --gpus all --mount type=bind,source=/storage,target=/app/storage";
#set AF2_CMD = "docker run --user "${u}":"${g}" --rm --gpus all --mount type=bind,source=/storage,target=/app/storage";
set prefix = "/app";
  endif

set n_preds = $4;

  if( -e $outdir2/${fasta:t:r}/ranking_debug.json && ! -z $outdir2/${fasta:t:r}/ranking_debug.json ) then
echo "Final file ranking_debug.json found. Aborting..." | tee -a $outdir.log;
echo "This probably means that this prediction was already calculated before." | tee -a $outdir.log;
echo "Moving $outdir2 to $outdir" | tee -a $outdir.log;
rm -f $outdir2/running.txt;
if($status) goto failed;
echo 0 >! $outdir2/finished.txt;
if($status) goto failed;
chown ${u}:${g} $outdir2/finished.txt
if($status) goto failed;
mv $outdir2 $outdir;
if($status) goto failed;
exit 0;
  endif

#the PROCINFO below set the array traversing order according to https://www.gnu.org/software/gawk/manual/html_node/Controlling-Scanning.html
awk 'function ltrim(s){sub(/^[ \t\r\n]+/,"",s);return s}function rtrim(s){sub(/[ \t\r\n]+$/,"",s);return s}function trim(s){return rtrim(ltrim(s))}BEGIN{i=0;seq=""}{if(substr($1,1,1)==">"){if(seq!=""){aarray[seq]++;if(aarray[seq]==1){aarray2[seq]=i;i++}};seq=""} else {seq=seq""toupper(trim($0))}}END{aarray[seq]++;if(aarray[seq]==1){aarray2[seq]=i;i++};PROCINFO["sorted_in"]="@val_num_asc";for(seq in aarray2){print seq;}}' $fasta >! $$.unique_seqs;
if($status) goto failed;
set nus = `cat $$.unique_seqs | wc -l`;
if($status) goto failed;

set seed = "";
  if($#argv > 4) then
    if( "$5" != "") then
      if(`echo "$5" | awk '{if(math($1,/[^0-9]+/)){print 1} else {print 0}}'`) then
echo "${caller}: ERROR 3rd argument is not a numeric seed" | tee -a $outdir.log;
goto failed;
      else
set seed = "--random_seed $5";
      endif
    endif
  endif

set skipqueryDBs = 0;
  if($#argv > 5) then
    if( "$6" != "") then
      if("$6" == "--skipqueryDBs") then
set skipqueryDBs = 1;
echo "Skipping query DBs" | tee -a $outdir.log;
      else
echo "${caller}: ERROR: 4th argument is not --skipqueryDBs" | tee -a $outdir.log;
goto failed;
      endif
    endif
  endif

set ignoretags = 0;
  if($#argv > 6) then
    if( "$7" != "") then
      if("$7" == "--ignoretags") then
set ignoretags = 1;
      else
echo "${caller}: ERROR: 5th argument is not --ignoretags" | tee -a $outdir.log;
goto failed;
      endif
    endif
  endif

#all possible templates
set templatedate = 3000-01-01
  if($#argv > 7) then
    if( "$8" != "") then
      if("$8" == "--notemplates") then
#no templates
set templatedate = 0001-01-01
      else
echo "${caller}: ERROR: 6th argument is not --notemplates" | tee -a $outdir.log;
goto failed;
      endif
    endif
  endif

set model = "";
  if($#argv > 8) then
    if( "$9" != "") then
      if( "$9" != "1" && "$9" != "2" && "$9" != "3" && "$9" != "4" && "$9" != "5" ) then
echo "${caller}: ERROR 7th argument is not a model NUMBER [1-5]" | tee -a $outdir.log;
goto failed;
      else
set model = "$9";
      endif
    endif
  endif

#relax is not working with the current setup (Debian 12 + nvidia drivers + cuda version)
#disabled for now
set common_args = "$seed --models_to_relax none --num_multimer_predictions_per_model $n_preds"
set data = "$prefix/storage/Alphafold/data";

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

  if( ! $ignoretags ) then
#Search for tags 
set i = 1;
set tagname = ( his-tag HQ-tag    HN-tag        HAT-tag           SpyTag        KTag      SnoopTag     SpyTag002     SnoopTagJr          DogTag             AVI-atg );
    foreach tag ( HHHHHH HQHQHQ HNHNHNHNHNHN KDHLIHNVHKEEHAHAHNK AHIVMVDAYKPTK ATHIKFSKRD KLGDIEFIKVNK VPTIVMVDAYKRYK KLGSIEFIKVNK DIPATYEFTDGKHYITNEPIPPK GLNDIFEAQKIEWHE )
grep -n $tag $$.unique_seqs;
      if( $status == 0 ) then
echo "The sequence(s) above seem(s) to have a $tagname[$i] $tag" | tee -a $outdir.log;
echo "Please verify that is what you want and run with the flag --ignoretags" | tee -a $outdir.log;
goto failed;
      endif
@ i++;
    end
  endif

if ($skipqueryDBs) goto labelskipqueryDBs;

set SeqWithoutMatch = 0;
set i = 1;
set ids = "";
set nhits = 0;
  while ($i <= $nus)
set cur_seq = `head -$i $$.unique_seqs | tail -1`;
#the commands below should be able to verify if a sequence is already present in the PDB
echo "Querying PDB server for EXPERIMENTAL sequence ${cur_seq}" | tee -a $outdir.log;
curl -s X GET --header 'Accept:application/json' 'https://search.rcsb.org/rcsbsearch/v2/query?json=%7B%0A%20%20%22query%22%3A%20%7B%0A%20%20%20%20%22type%22%3A%20%22group%22%2C%0A%20%20%20%20%22logical_operator%22%3A%20%22and%22%2C%0A%20%20%20%20%22nodes%22%3A%20%5B%20%7B%0A%20%20%20%20%20%20%22type%22%3A%22terminal%22%2C%0A%20%20%20%20%20%20%22service%22%3A%22seqmotif%22%2C%0A%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%22value%22%3A%22'${cur_seq}'%22%2C%0A%20%20%20%20%20%20%20%20%22pattern_type%22%3A%22simple%22%2C%0A%20%20%20%20%20%20%20%20%22sequence_type%22%3A%22protein%22%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%7D%20%5D%0A%20%20%7D%2C%0A%20%20%22return_type%22%3A%20%22polymer_entity%22%2C%0A%20%20%22request_options%22%3A%20%7B%0A%20%20%20%20%22results_content_type%22%3A%20%5B%22experimental%22%5D%2C%0A%20%20%20%20%22scoring_strategy%22%3A%20%22combined%22%0A%20%20%7D%0A%7D' | awk -v FS=',' '{for(i=1; i<=NF; i++) print $i}' >! $$.fetchPDB.$i;
    if( ! -z $$.fetchPDB.$i ) then
grep 'identifier' $$.fetchPDB.$i | sed -e 's/"result_set":\[//' >! $$.fetchPDB_ids.$i;
grep 'score' $$.fetchPDB.$i >! $$.fetchPDB_scores.$i;
      if( ! -z $$.fetchPDB_ids.$i && ! -z $$.fetchPDB_scores.$i ) then
paste -d ' ' $$.fetchPDB_ids.$i $$.fetchPDB_scores.$i | sed -e 's/:/ /g' >! $$.fetchPDB.$i;
        if( `sort -n -k 4 $$.fetchPDB.$i | tail -1 | awk '{if($4>=0.999){print 1} else {print 0}}'` ) then
@ nhits++
echo "Sequence $i match with PDB https://www.rcsb.org/structure/"`sort -n -k 4 $$.fetchPDB.$i | tail -1 | awk '{print substr($2,2,4)}'` | tee -a $outdir.log;
echo -n " $i" >>! $$.matchPDB.`sort -n -k 4 $$.fetchPDB.$i | tail -1 | awk '{print substr($2,2,4)}'`;
set ids = "$ids "$i;
        else
echo "Sequence $i did not find high-score matches" | tee -a $outdir.log;
        endif
      else
echo "Matches for sequence $i returned no identifier/score values" | tee -a $outdir.log;
      endif
rm -f $$.fetchPDB.$i $$.fetchPDB_ids.$i $$.fetchPDB_scores.$i;
    else
echo "Sequence $i did not find any match" | tee -a $outdir.log;
set SeqWithoutMatch = 1;
rm -f $$.fetchPDB.$i;
    endif
@ i++
  end
  if( ! $SeqWithoutMatch ) then
    if( "$ids" != "" ) then
      foreach f ( $$.matchPDB.???? )
echo "" >> $f;
        if ( "`cat $f`" == "$ids" ) then
echo "All the requested sequences match the same PDB: https://www.rcsb.org/structure/${f:e} Aborting..." | tee -a $outdir.log;
rm -f $$.matchPDB.????;
goto failed;
        endif
      end
    endif
  endif

rm -f $$.matchPDB.????;

set SeqWithoutMatch = 0;
set i = 1;
set ids = "";
set nhits = 0;
  while ($i <= $nus)
set cur_seq = `head -$i $$.unique_seqs | tail -1`;
#the commands below should be able to verify if a sequence is already present in the AFDB
echo "Querying PDB server for COMPUTATIONAL sequence ${cur_seq}" | tee -a $outdir.log;
curl -s X GET --header 'Accept:application/json' 'https://search.rcsb.org/rcsbsearch/v2/query?json=%7B%0A%20%20%22query%22%3A%20%7B%0A%20%20%20%20%22type%22%3A%20%22group%22%2C%0A%20%20%20%20%22logical_operator%22%3A%20%22and%22%2C%0A%20%20%20%20%22nodes%22%3A%20%5B%20%7B%0A%20%20%20%20%20%20%22type%22%3A%22terminal%22%2C%0A%20%20%20%20%20%20%22service%22%3A%22seqmotif%22%2C%0A%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%22value%22%3A%22'${cur_seq}'%22%2C%0A%20%20%20%20%20%20%20%20%22pattern_type%22%3A%22simple%22%2C%0A%20%20%20%20%20%20%20%20%22sequence_type%22%3A%22protein%22%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%7D%20%5D%0A%20%20%7D%2C%0A%20%20%22return_type%22%3A%20%22polymer_entity%22%2C%0A%20%20%22request_options%22%3A%20%7B%0A%20%20%20%20%22results_content_type%22%3A%20%5B%22computational%22%5D%2C%0A%20%20%20%20%22scoring_strategy%22%3A%20%22combined%22%0A%20%20%7D%0A%7D' | awk -v FS=',' '{for(i=1; i<=NF; i++) print $i}' >! $$.fetchAFDB.$i;
    if( ! -z $$.fetchAFDB.$i ) then
grep 'identifier' $$.fetchAFDB.$i | sed -e 's/"result_set":\[//' >! $$.fetchAFDB_ids.$i;
grep 'score' $$.fetchAFDB.$i >! $$.fetchAFDB_scores.$i;
      if( ! -z $$.fetchAFDB_ids.$i && ! -z $$.fetchAFDB_scores.$i ) then
paste -d ' ' $$.fetchAFDB_ids.$i $$.fetchAFDB_scores.$i | sed -e 's/:/ /g' >! $$.fetchAFDB.$i;
        if( `sort -n -k 4 $$.fetchAFDB.$i | tail -1 | awk '{if($4>=0.999){print 1} else {print 0}}'` ) then
@ nhits++
echo "Sequence $i match with AFDB https://www.rcsb.org/structure/"`sort -n -k 4 $$.fetchAFDB.$i | tail -1 | awk '{print substr($2,2,13)}'` | tee -a $outdir.log;
echo -n " $i" >>! $$.matchAFDB.`sort -n -k 4 $$.fetchAFDB.$i | tail -1 | awk '{print substr($2,2,13)}'`;
set ids = "$ids "$i;
        else
echo "Sequence $i did not find high-score matches" | tee -a $outdir.log;
        endif
      else
echo "Matches for sequence $i returned no identifier/score values" | tee -a $outdir.log;
      endif
rm -f $$.fetchAFDB.$i $$.fetchAFDB_ids.$i $$.fetchAFDB_scores.$i;
    else
echo "Sequence $i did not find any match" | tee -a $outdir.log;
set SeqWithoutMatch = 1;
rm -f $$.fetchAFDB.$i;
    endif
@ i++
  end
  if( ! $SeqWithoutMatch ) then
    if( "$ids" != "" ) then
      foreach f ( $$.matchAFDB.?????????? )
echo "" >> $f;
        if ( "`cat $f`" == "$ids" ) then
echo "All the requested sequences match the same AFDB: https://alphafold.ebi.ac.uk/entry/${f:e} Aborting..." | tee -a $outdir.log;
rm -f $$.matchAFDB.??????????;
goto failed;
        endif
      end
    endif
  endif

rm -f $$.matchAFDB.??????????;

labelskipqueryDBs:

set alphabet = ( A B C D E F G H I J K L M N O P Q R S T U V W X Y Z );

#AF2seqs.txt has 1st column: sequence 2nd column: msa_folder/ with the .sto/.a3m files inside.
  if (! -e /storage/Data/AF2jobs/AF2seqs.txt) then
touch /storage/Data/AF2jobs/AF2seqs.txt;
if($status) goto failed;
chmod 664 /storage/Data/AF2jobs/AF2seqs.txt;
if($status) goto failed;
  else
set i = 1;
if($status) goto failed;
    foreach seq ( `cat $$.unique_seqs` )
      if( `awk -v seq=$seq 'BEGIN{p=0}{if(seq==$1) {print 1; p=1; exit;}}END{if(p==0) {print 0} }' /storage/Data/AF2jobs/AF2seqs.txt` ) then
mkdir -m 775 -p $outdir2/${fasta:t:r}/msas/$alphabet[$i]/;
if($status) goto failed;
chmod 775 $outdir2/${fasta:t:r}/msas/;
if($status) goto failed;
chmod 775 $outdir2/${fasta:t:r}/;
if($status) goto failed;
chmod 775 $outdir2/;
if($status) goto failed;
chown -R ${u}:${g} $outdir2/;
if($status) goto failed;
cp -a `awk -v seq=$seq '{if(seq==$1) {print $2"/*"}}' /storage/Data/AF2jobs/AF2seqs.txt` $outdir2/${fasta:t:r}/msas/$alphabet[$i]/
      endif
@ i++
    end
  endif

cp -a $fasta /storage/Data/AF2jobs/${fasta:t};
if($status) goto failed;

set i = 1;
set nhits = 0;
echo "Job length: $len; #UniqueSequences: $nus" | tee -a $outdir.log;
#optimized for parallel models by LFE (stopat 0 use no GPU):
  if($RUN_PARALLEL) then
setenv XLA_PYTHON_CLIENT_PREALLOCATE false
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 1.0
set stopat = "--stopat 0";
  else
setenv XLA_PYTHON_CLIENT_PREALLOCATE true
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 100.0
set stopat = "";
  endif
  if( $INSTALLATION_TYPE == "docker" ) then
eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=1.0 alphafold --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas $stopat" |& tee -a $outdir.log;
if($status) goto failed;
  else
$AF2_CMD --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas $stopat |& tee -a $outdir.log;
if($status) goto failed;
  endif

  if(! $RUN_PARALLEL) then
goto legrandfinale;
  endif

  if( $model == "") then
    if($len < 500) then
#try 5 processes in parallel
#optimized for parallel models by LFE:
setenv XLA_PYTHON_CLIENT_PREALLOCATE false
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 20.0
set i = 1;
      while ( $i < 6 )
        if( $INSTALLATION_TYPE == "docker" ) then
( eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=20.0 alphafold --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" & ) |& tee -a $outdir.$i.log &;
        else
( $AF2_CMD --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i & ) |& tee -a $outdir.$i.log &;
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
( eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=33.333 alphafold --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" & ) |& tee -a $outdir.$i.log &;
         else
( $AF2_CMD --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i & ) |& tee -a $outdir.$i.log &;
         endif
@ i++
        end
wait;
#try 2 processes in parallel
setenv XLA_PYTHON_CLIENT_PREALLOCATE false
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 50.0
        while ( $i < 6 )
          if( $INSTALLATION_TYPE == "docker" ) then
( eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=50.0 alphafold --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" & ) |& tee -a $outdir.$i.log &;
          else
( $AF2_CMD --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i & ) |& tee -a $outdir.$i.log &;
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
( eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=50.0 alphafold --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" & ) |& tee -a $outdir.$i.log &;
            else
( $AF2_CMD --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i & ) |& tee -a $outdir.$i.log &;
            endif
@ i++
          end
wait;
          while ( $i < 5 )
            if( $INSTALLATION_TYPE == "docker" ) then
( eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env XLA_PYTHON_CLIENT_MEM_FRACTION=50.0 alphafold --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" & ) |& tee -a $outdir.$i.log &;
            else
( $AF2_CMD --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i & ) |& tee -a $outdir.$i.log &;
            endif
@ i++
          end
wait;
#agressive memory use
setenv XLA_PYTHON_CLIENT_PREALLOCATE true
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 100.0
          if( $INSTALLATION_TYPE == "docker" ) then
eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=true --env XLA_PYTHON_CLIENT_MEM_FRACTION=100.0 alphafold --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i" |& tee -a $outdir.$i.log;
          else
$AF2_CMD --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $i |& tee -a $outdir.$i.log;
          endif
        endif
      endif
    endif
  else
#agressive memory use
setenv XLA_PYTHON_CLIENT_PREALLOCATE true
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 100.0
    if( $INSTALLATION_TYPE == "docker" ) then
eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=true --env XLA_PYTHON_CLIENT_MEM_FRACTION=100.0 alphafold --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $model" |& tee -a $outdir.log;
    else
$AF2_CMD --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $model |& tee -a $outdir.log;
    endif
  endif
  if($len < 500) then
set i = 1;
    while ( $i < 6 )
cat $outdir.$i.log >>! $outdir.log;
rm -f $outdir.$i.log;
@ i++;
    end
  endif
#unexpected good thing: if any of the run-in-background inferences above fail because lack of vRAM
#the run_alphafold below run the missing ones one-by-one
#agressive memory use
setenv XLA_PYTHON_CLIENT_PREALLOCATE true
setenv XLA_PYTHON_CLIENT_MEM_FRACTION 100.0

set final_stopat = 6;
  if($model != "") then
set final_stopat = 7;
  endif

  if( $INSTALLATION_TYPE == "docker" ) then
eval $AF2_CMD" --env XLA_PYTHON_CLIENT_PREALLOCATE=true --env XLA_PYTHON_CLIENT_MEM_FRACTION=100.0 alphafold --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $final_stopat" |& tee -a $outdir.log;
  else
$AF2_CMD --fasta_paths $prefix/storage/Data/AF2jobs/${fasta:t} --output_dir ${prefix}$outdir2 --data_dir $data --db_preset=$msa_dbs_mode --uniref90_database_path $data/uniref90/uniref90.fasta --mgnify_database_path $data/mgnify/mgy_clusters.fa --template_mmcif_dir $data/pdb_mmcif/mmcif_files --max_template_date=$templatedate --obsolete_pdbs_path $data/pdb_mmcif/obsolete.dat --use_gpu_relax=True --model_preset=multimer $preset_dependent_args $msa_dbs_mode_dependent_args --use_precomputed_msas --stopat $final_stopat |& tee -a $outdir.log;
  endif
if($status) goto failed;

legrandfinale:

chown ${u}:${g} $outdir.log;
if($status) goto failed;

set i = 1;
if($status) goto failed;
  foreach seq ( `cat $$.unique_seqs` )
grep -i "$seq" /storage/Data/AF2jobs/AF2seqs.txt >& /dev/null;
    if($status == 1) then #match not found
set tmpdir = `mktemp -u -d -t -p /storage/Data/AF2jobs/AF2msas/`;
if($status) goto failed;
mkdir -m 775 -p $tmpdir;
if($status) goto failed;
chown -R ${u}:${g} $tmpdir/;
if($status) goto failed;
echo "$seq $tmpdir" >>! /storage/Data/AF2jobs/AF2seqs.txt;
if($status) goto failed;
cp -a $outdir2/${fasta:t:r}/msas/$alphabet[$i]/* $tmpdir/; 
if($status) goto failed;
    endif
@ i++
  end

  if(-e /storage/Data/AF2jobs/AF2.log ) then
grep $jobname /storage/Data/AF2jobs/AF2.log | awk -v FS="<td>" '{print $3}' | sed -e 's/<\/td>//' >& /dev/null;
    if($status) then
echo "${caller}: grep jobname failed." |& tee -a $outdir.log;
    else
set email = `grep $jobname /storage/Data/AF2jobs/AF2.log | awk -v FS="<td>" '{print $3}' | sed -e 's/<\/td>//'`
set hostname = `hostname`;
set email = "${email}@ibs.fr"
sudo -u \#50809 mail -s '[no-reply] AF2 results' -c leandro.estrozi@ibs.fr -- $email << EOF
AF2 results are ready.
You can downdload then at: http://${hostname}.ibs.fr:8080/AF2IBS/browse/${jobname}?filepath=
EOF
if($status) echo "${caller}: send mail failed." |& tee -a $outdir.log;
    endif
  endif

rm -f $outdir2/running.txt;
if($status) goto failed;
echo 0 >! $outdir2/finished.txt;
if($status) goto failed;
chown ${u}:${g} $outdir2/finished.txt
if($status) goto failed;
echo "Moving $outdir2 to $outdir" |& tee -a $outdir.log;
mv $outdir2 $outdir;
if($status) goto failed;
chown -R ${u}:${g} $outdir;
if($status) goto failed;
rm -f /storage/Data/AF2jobs/${fasta:t};
if($status) goto failed;
rm -f $$.unique_seqs;
if($status) goto failed;

exit 0;

failed:
rm -f $outdir2/running.txt;
echo "last cmd before fail: $_" >! $outdir2/failed.txt;
chown ${u}:${g} $outdir2/failed.txt
chown ${u}:${g} $outdir.log
echo "Moving $outdir2 to $outdir" |& tee -a $outdir.log;
mv $outdir2 $outdir;
rm -f /storage/Data/AF2jobs/${fasta:t};
rm -f $$.unique_seqs;
exit 1;
