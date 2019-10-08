
use List::Util qw( min max );
our %AA3TO1 = qw(ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V);
our %AA1TO3 = reverse %AA3TO1;

sub permute (&@);  # predeclare sub, paste sub at bottom
 
if (@ARGV <4 )
{
	die "Error: need seven parameters: domain_list, domain model folder, query file(fasta), target id, output dir, modeller path, model number. \n";
}
#### 
$input_pdb = shift @ARGV;
$fasta_file = shift @ARGV;
$outputfolder = shift @ARGV;
$proc_num = shift @ARGV;
$fragments_prior = shift @ARGV;

-d $outputfolder || `mkdir $outputfolder`;

$fragseq = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
#get query name and sequence 
open(FASTA, $fasta_file) || die "Error: can't read fasta file $fasta_file.\n";
@content = <FASTA>;
$target_id = shift @content;
chomp $target_id;
$qseq = shift @content;
chomp $qseq;
close FASTA;



#rewrite fasta file if it contains lower-case letter
if ($qseq =~ /[a-z]/)
{
	print "There are lower case letters in the input file. Convert them to upper case.\n";
	$qseq = uc($qseq);
	open(FASTA, ">$outputfolder/$target_id.fasta") || die "Error: can't rewrite fasta file.\n";
	print FASTA "$target_id\n$qseq\n";
	close FASTA;
}


if ($target_id =~ /^>/)
{
	$target_id = substr($target_id, 1); 
}
else
{
	die "Error: fasta foramt error.\n"; 
}

@prior_fragment_orders = ();
if(-e $fragments_prior)
{
  
  open(TMP, $fragments_prior) || die "Error: can't read fasta file $fragments_prior.\n";
  @content = <TMP>;
  foreach $line (@content)
  {
    @tmp = split(/\s+/,$line);
    if(@tmp<2)
    {
      next;
    }
    $new = join(' ',@tmp);
    push @prior_fragment_orders,$new;
    print "Loading prior fragment orders: $new!\n";
    @tmp_new = reverse @tmp;
    $new = join(' ',@tmp_new);
    push @prior_fragment_orders,$new;
    print "Loading prior fragment orders: $new!\n";
  }
  close TMP;
}

#####  (1) initialize sequence in fragment


open INPUTPDB, $input_pdb or die "ERROR! Could not open $input_pdb";
@lines_PDB = <INPUTPDB>;
close INPUTPDB;


-d "$outputfolder/frag_dir" || `mkdir $outputfolder/frag_dir`;

@PDB_temp=();
$frag_num=0;
$frag_start = 0;
%frag_CAs = ();
$CA_num = 0;
foreach (@lines_PDB) {
	$line = $_;
	chomp $line;
 
  if(substr($line,0,3) eq 'TER')
  {
    $frag_start = 0;
    next;  
  }
  
  
	next if $line !~ m/^ATOM/;
	@tmp = split(/\s+/,$line);
  $atomtype = parse_pdb_row($line,"aname");
  next if $atomtype != 'CA';
  
  
  if($frag_start == 0)
  {
    $frag_num++;
    $frag_start = 1;
    if($frag_num>1)
    {
      $idx = $frag_num-1;
      $fragindx = substr($fragseq,$idx-1,1);
      @keys_f = keys %frag_CAs;
      $frag_CAs{$fragindx} = $CA_num;
      @keys_f = keys %frag_CAs;
      close TMP;
    }
    
    $fragindx = substr($fragseq,$frag_num-1,1);
    $fragfile = "$outputfolder/frag_dir/frag${fragindx}_tmp.pdb";
    open(TMP,">$fragfile") || die "Failed to open file $fragfile\n\n";
    print TMP "$line\n";
    $CA_num = 1;
    next;
    
  }
  
  print TMP "$line\n";
  $CA_num ++;
  
}
$fragindx = substr($fragseq,$frag_num-1,1);
$frag_CAs{$fragindx} = $CA_num;
close TMP;

print "Total $frag_num fragments are found\n";

### sort fragments by length
$idx=0;
foreach $fragidx (sort { $frag_CAs{$b} <=> $frag_CAs{$a} } keys %frag_CAs) 
{
  $frag_len = $frag_CAs{$fragidx};
  if($frag_len<3)
  {
    $frag_num = $frag_num-1;
    next;
  }
  $frag1 = "$outputfolder/frag_dir/frag${fragidx}_tmp.pdb";

  ##### (3) add pulchar and side-chain
  
  open(TMP1, "$frag1") || die "Failed to open1 $frag1\n";
  @content1 = <TMP1>;
  close TMP1;
  $Ca_index = 0;
  $atom_index=0;  
  $this_rchain="";
  $pdb_seq="";
  $idx++;
  $fragidx_new = substr($fragseq,$idx-1,1);
  $frag1_sort = "$outputfolder/frag_dir/frag${fragidx_new}.pdb";
  open(OUTPDB,">$frag1_sort") || die "Failed to open $frag1_sort\n";
  foreach(@content1)
  {
    	$line = $_;
    	chomp $line;
    	next if $line !~ m/^ATOM/;
    	$atomCounter = parse_pdb_row($line,"anum");
    	$atomtype = parse_pdb_row($line,"aname");
    	$resname = parse_pdb_row($line,"rname");
    	$chainid = parse_pdb_row($line,"chain");
    	$resCounter = parse_pdb_row($line,"rnum");
      $this_rchain = $chainid;
      
      if($atomtype eq 'CA')
      {
        $Ca_index++;
        $pdb_seq .= $AA3TO1{$resname};
      }
      $atom_index++;
       
      $x = parse_pdb_row($line,"x");
      $y = parse_pdb_row($line,"y");
      $z = parse_pdb_row($line,"z");
     
    	
      
    	$rnum_string = sprintf("%4s", $Ca_index);
    	$anum_string = sprintf("%5s", $atom_index);
    	$atomtype = sprintf("%4s", $atomtype);
    	$x = sprintf("%8s", $x);
    	$y = sprintf("%8s", $y);
    	$z = sprintf("%8s", $z);
    	$row = "ATOM  ".$anum_string.$atomtype."  ".$resname." ".$chainid.$rnum_string."    ".$x.$y.$z."\n";
    	print OUTPDB $row;
  }
  close OUTPDB;
  `rm $frag1`;
  print "Fragement $fragidx_new has $frag_len atoms\n";
}

### get forward and reverse pdb

for($idx=1;$idx<=$frag_num;$idx++)
{
  $idx1 = $idx;
  $fragidx = substr($fragseq,$idx-1,1);
  $frag1 = "$outputfolder/frag_dir/frag$fragidx.pdb";
  $frag1_f = "$outputfolder/frag_dir/frag${fragidx}_f.pdb";
  $frag1_r = "$outputfolder/frag_dir/frag${fragidx}_r.pdb";
  open OUTPDB, ">$frag1_f" or die "ERROR! Could not open $frag1_f\n";
  
  open INPUTPDB, "$frag1" or die "ERROR! Could not open $frag1";
  my @lines_PDB = <INPUTPDB>;
  close INPUTPDB;
  
  @PDB_temp=();
  foreach (@lines_PDB) {
  	$line = $_;
  	chomp $line;
  	next if $line !~ m/^ATOM/;
  
  	$atomCounter = parse_pdb_row($line,"anum");
  	$atomtype = parse_pdb_row($line,"aname");
  	$resname = parse_pdb_row($line,"rname");
  	$chainid = parse_pdb_row($line,"chain");
  	$resCounter = parse_pdb_row($line,"rnum");
    $this_rchain = $chainid;
   
    if($atomtype eq 'CA')
    {
      push @PDB_temp,$line;
    }
    
    if($resCounter<=length($qseq))
    {
      $newresname = $AA1TO3{substr($qseq,$resCounter-1,1)};
    }
    
    $x = parse_pdb_row($line,"x");
    $y = parse_pdb_row($line,"y");
    $z = parse_pdb_row($line,"z");
  
  	my $rnum_string = sprintf("%4s", $resCounter);
  	my $anum_string = sprintf("%5s", $atomCounter);
  	my $atomtype = sprintf("%4s", $atomtype);
  	my $x = sprintf("%8s", $x);
  	my $y = sprintf("%8s", $y);
  	my $z = sprintf("%8s", $z);
  	my $row = "ATOM  ".$anum_string.$atomtype."  ".$newresname." ".$chainid.$rnum_string."    ".$x.$y.$z."\n";
  	print OUTPDB $row;
  }
  print OUTPDB "END\n";
  close OUTPDB;
  
  
  
  
  #####  (2) initialize sequence in reverse fragment
  
  open OUTPDB, ">$frag1_r" or die "ERROR! Could not open $frag1_r\n";
  
  @PDB_temp = reverse(@PDB_temp);
  $Ca_index = 0;
  $atom_index=0;
  foreach (@PDB_temp) {
  	$line = $_;
  	chomp $line;
  	next if $line !~ m/^ATOM/;
  	$atomCounter = parse_pdb_row($line,"anum");
  	$atomtype = parse_pdb_row($line,"aname");
  	$resname = parse_pdb_row($line,"rname");
  	$chainid = parse_pdb_row($line,"chain");
  	$resCounter = parse_pdb_row($line,"rnum");
    $this_rchain = $chainid;
    next if $atomtype != 'CA';
    
    $Ca_index++;
    $atom_index++;
     
   
    $x = parse_pdb_row($line,"x");
    $y = parse_pdb_row($line,"y");
    $z = parse_pdb_row($line,"z");
  	
  	my $rnum_string = sprintf("%4s", $Ca_index);
  	my $anum_string = sprintf("%5s", $atom_index);
  	my $atomtype = sprintf("%4s", $atomtype);
  	my $x = sprintf("%8s", $x);
  	my $y = sprintf("%8s", $y);
  	my $z = sprintf("%8s", $z);
  	my $row = "ATOM  ".$anum_string.$atomtype."  ".$resname." ".$chainid.$rnum_string."    ".$x.$y.$z."\n";
  	print OUTPDB $row;
  }
  print OUTPDB "END\n";
  close OUTPDB; 
}


#### generate combination of fragments
@frag_list = ();
for($idx=1;$idx<=$frag_num;$idx++)
{
  $fragindx = substr($fragseq,$idx-1,1);
  push @frag_list,"frag$fragindx";
}

$count=0;

if(-d "$outputfolder/frag_merge_dir/")
{
  `rm -rf $outputfolder/frag_merge_dir/*`;
}else{
  `mkdir -p $outputfolder/frag_merge_dir/`;
}

@permuation_frags = ();
permute{ 
  $new_list = join(' ',@_);
  push @permuation_frags, $new_list;
  
  
} @frag_list;

%frag_directions = ();
@running_files = ();
foreach (@permuation_frags)
{
  if(@prior_fragment_orders>0)
  {
      $found = 0;
      foreach $item (@prior_fragment_orders)
      {
        if(index($_,"$item")>=0)
        {
          $found = 1;
        }
      }
      if($found != 1)
      {
        next;
      }
  }
  
  @frag_orders = split(' ',$_);
  ############## need merge fragments, check possible amino acids between fragments
  #@frag_orders = qw(frag8 frag10 frag5 frag6 frag1 frag2 frag4 frag3 frag9 frag7);
  %added_gap=();
  $large_gaps=0;
  for($idx=0;$idx<$frag_num-1;$idx++)
  {
    $idx1 = $idx;
    $idx2 = $idx+1;
    $frag1 = "$outputfolder/frag_dir/".$frag_orders[$idx].".pdb";
    $frag2 = "$outputfolder/frag_dir/".$frag_orders[$idx+1].".pdb";
    
    open(TMP1, "$frag1") || die "Failed to open $frag1\n";
    @content1 = <TMP1>;
    close TMP1;
    
    open(TMP2, "$frag2") || die "Failed to open $frag2\n";
    @content2 = <TMP2>;
    close TMP2;
    
    $line1_first = shift @content1;
    $line1_last = pop @content1;
    $line2_first = shift @content2;
    $line2_last = pop @content2;
    
    $line1_first_x = parse_pdb_row($line1_first,"x");
    $line1_first_y = parse_pdb_row($line1_first,"y");
    $line1_first_z = parse_pdb_row($line1_first,"z");
    
    $line1_last_x = parse_pdb_row($line1_last,"x");
    $line1_last_y = parse_pdb_row($line1_last,"y");
    $line1_last_z = parse_pdb_row($line1_last,"z");
    
    $line2_first_x = parse_pdb_row($line2_first,"x");
    $line2_first_y = parse_pdb_row($line2_first,"y");
    $line2_first_z = parse_pdb_row($line2_first,"z");
    
    $line2_last_x = parse_pdb_row($line2_last,"x");
    $line2_last_y = parse_pdb_row($line2_last,"y"); 
    $line2_last_z = parse_pdb_row($line2_last,"z");
    
    
    $distance1 = sqrt(($line1_first_x-$line2_first_x)*($line1_first_x-$line2_first_x) + ($line1_first_y-$line2_first_y)*($line1_first_y-$line2_first_y)+($line1_first_z-$line2_first_z)*($line1_first_z-$line2_first_z));
     
    $distance2 = sqrt(($line1_last_x-$line2_first_x)*($line1_last_x-$line2_first_x) + ($line1_last_y-$line2_first_y)*($line1_last_y-$line2_first_y)+($line1_last_z-$line2_first_z)*($line1_last_z-$line2_first_z));
    
    $distance3 = sqrt(($line1_first_x-$line2_last_x)*($line1_first_x-$line2_last_x) + ($line1_first_y-$line2_last_y)*($line1_first_y-$line2_last_y)+($line1_first_z-$line2_last_z)*($line1_first_z-$line2_last_z));
    
    $distance4 = sqrt(($line1_last_x-$line2_last_x)*($line1_last_x-$line2_last_x) + ($line1_last_y-$line2_last_y)*($line1_last_y-$line2_last_y)+($line1_last_z-$line2_last_z)*($line1_last_z-$line2_last_z));
    

    @tmp =();
    
  
    if(exists($frag_directions{$frag_orders[$idx]}) and $frag_directions{$frag_orders[$idx]} eq 'Reverse')
    {
      push @tmp,$distance1;
      push @tmp,$distance3;
    }elsif(exists($frag_directions{$frag_orders[$idx]}) and $frag_directions{$frag_orders[$idx]} eq 'Forward')
    {
      push @tmp,$distance2;
      push @tmp,$distance4;
    }else{
      push @tmp,$distance1;
      push @tmp,$distance2;
      push @tmp,$distance3;
      push @tmp,$distance4;
    }
    @dist_array = sort { $a <=> $b } @tmp;
    
    $distance = $dist_array[0];
    
    if($distance == $distance1)
    {
      print "$frag_orders[$idx].pdb is reverse and $frag_orders[$idx+1].pdb is forward\n";
      if(exists($frag_directions{$frag_orders[$idx]}) and $frag_directions{$frag_orders[$idx]} ne 'Reverse') ### means conflict with previous one, need pass
      {
        $large_gaps = 1000;
      }
      $frag_directions{$frag_orders[$idx]} = 'Reverse';
      $frag_directions{$frag_orders[$idx+1]} = 'Forward';
    }elsif($distance == $distance2)
    {
      print "$frag_orders[$idx].pdb is forward and $frag_orders[$idx+1].pdb is forward\n";
      
      if(exists($frag_directions{$frag_orders[$idx]}) and $frag_directions{$frag_orders[$idx]} ne 'Forward') ### means conflict with previous one, need pass
      {
        $large_gaps = 1000;
      }
      $frag_directions{$frag_orders[$idx]} = 'Forward';
      $frag_directions{$frag_orders[$idx+1]} = 'Forward';
    }elsif($distance == $distance3)
    {
      print "$frag_orders[$idx].pdb is reverse and $frag_orders[$idx+1].pdb is reverse\n";
      
      if(exists($frag_directions{$frag_orders[$idx]}) and $frag_directions{$frag_orders[$idx]} ne 'Reverse') ### means conflict with previous one, need pass
      {
        $large_gaps = 1000;
      }
      $frag_directions{$frag_orders[$idx]} = 'Reverse';
      $frag_directions{$frag_orders[$idx+1]} = 'Reverse';
    }elsif($distance == $distance4)
    {
      print "$frag_orders[$idx].pdb is forward and $frag_orders[$idx+1].pdb is reverse\n";
      
      if(exists($frag_directions{$frag_orders[$idx]}) and $frag_directions{$frag_orders[$idx]} ne 'Forward') ### means conflict with previous one, need pass
      {
        $large_gaps = 1000;
      }
      $frag_directions{$frag_orders[$idx]} = 'Forward';
      $frag_directions{$frag_orders[$idx+1]} = 'Reverse';
    }else{
      die "Wrong distance $distance and $tmp\n";
    }
    
    #print("Minimum is $distance\n");
    $num_aa = int($distance/3.8);
    if($num_aa > 5) ## assume missing at most 5 CA
    {
      $large_gaps ++;
    }
    #print "Distance: $distance Possible aa: $num_aa\n\n";
    
    $added_gap{$idx} = $num_aa;
  }
  if($large_gaps > $frag_num/2)
  {
    #print "Has large gap, Pass\n";
    next;
  }
  
  ## get CA-atom number
  open(OUTPDB,">$outputfolder/frag_dir/frag_merged.pdb") || die "Failed to open $outputfolder/frag_dir/\n";

  $Ca_index = 0;
  $atom_index=0;  
  
  $pdb_seq="";
  $pdb_seq_align="";
  $this_rchain="";
  for($idx=0;$idx<$frag_num;$idx++)
  {
    $idx1 = $idx;
    
    $frag1 = "$outputfolder/frag_dir/".$frag_orders[$idx]."_f.pdb";
    if($frag_directions{$frag_orders[$idx]} eq 'Reverse')
    {
      $frag1 = "$outputfolder/frag_dir/".$frag_orders[$idx]."_r.pdb";
    }
    open(TMP1, "$frag1") || die "Failed to open $frag1\n";
    @content1 = <TMP1>;
    close TMP1;
    
    foreach(@content1)
    {
      	$line = $_;
      	chomp $line;
      	next if $line !~ m/^ATOM/;
      	$atomCounter = parse_pdb_row($line,"anum");
      	$atomtype = parse_pdb_row($line,"aname");
      	$resname = parse_pdb_row($line,"rname");
      	$chainid = parse_pdb_row($line,"chain");
      	$resCounter = parse_pdb_row($line,"rnum");
        $this_rchain = $chainid;
        next if $atomtype != 'CA';
        
        $Ca_index++;
        $atom_index++;
         
        if($Ca_index<=length($qseq))
        {
          $newresname = $AA1TO3{substr($qseq,$Ca_index-1,1)};
        }
       
        $x = parse_pdb_row($line,"x");
        $y = parse_pdb_row($line,"y");
        $z = parse_pdb_row($line,"z");
       
      	
        $pdb_seq .= $AA3TO1{$newresname};
        $pdb_seq_align .= $AA3TO1{$newresname};
      	my $rnum_string = sprintf("%4s", $Ca_index);
      	my $anum_string = sprintf("%5s", $atom_index);
      	my $atomtype = sprintf("%4s", $atomtype);
      	my $x = sprintf("%8s", $x);
      	my $y = sprintf("%8s", $y);
      	my $z = sprintf("%8s", $z);
      	my $row = "ATOM  ".$anum_string.$atomtype."  ".$newresname." ".$chainid.$rnum_string."    ".$x.$y.$z."\n";
      	print OUTPDB $row;
    }
    if($idx<$frag_num)
    {
      $Ca_index = $Ca_index +$added_gap{$idx};
      $pdb_seq .= &generate_gaps($added_gap{$idx});
      $pdb_seq_align .= &generate_aa($added_gap{$idx});
      
    }
    
    
  }
  close OUTPDB;
  
  ## check if the total CA larger than original sequence
  if(length($pdb_seq_align)>length($qseq))
  {
    next;
  }
  #print("$pdb_seq_align\n$pdb_seq\n\n");
  #print("frag_merged_length: ".length($pdb_seq_align)."\n");
  #print("original seq length: ".length($qseq)."\n");
  $count++;
  

  `/scratch/jh7x3/multicom/tools/pulchra304/pulchra -c -s $outputfolder/frag_dir/frag_merged.pdb`;
  
  $init_pdb = "$outputfolder/frag_dir/frag_merged.pdb";
  if(!(-e "$outputfolder/frag_dir/frag_merged.rebuilt.pdb"))
  {
    $init_pdb = "$outputfolder/frag_dir/frag_merged.pdb";
  }else{
  
    `/scratch/jh7x3/multicom/tools/scwrl4/Scwrl4 -i $outputfolder/frag_dir/frag_merged.rebuilt.pdb -o  $outputfolder/frag_dir/frag_merged_scwrl.pdb`;
    if(!(-e "$outputfolder/frag_dir/frag_merged_scwrl.pdb"))
    {
      $init_pdb = "$outputfolder/frag_dir/frag_merged.rebuilt.pdb";
    }else{
      $init_pdb = "$outputfolder/frag_dir/frag_merged_scwrl.pdb";
    }
  }
  
    `mkdir -p $outputfolder/frag_merge_dir/frag_merged$count/`;
    `cp $init_pdb $outputfolder/frag_merge_dir/frag_merged$count/aaaaA.atm`;
    print "Found potential: @frag_orders\n";
    print "$pdb_seq\n$pdb_seq_align\n\n";
    
    #### start build models
    $pir_file = "$outputfolder/frag_merge_dir/frag_merged$count.pir";
    open(PIR, ">$pir_file") || die "can't create pir file $pir_file.\n";

    $dlen = length($pdb_seq);
    print PIR "C;cover size:X; local alignment length=X (original info = aaaaA	X	X	0	0)\n";
    print PIR ">P1;aaaaA\n";
    print PIR "structureX:aaaaA: 1: $this_rchain: $dlen: : : : : \n"; 
    print PIR "$pdb_seq*\n\n";
    
    print PIR "C;query_length:X num_of_temps:X cover_ratio:X cover:X not_cover:X\n"; 
    print PIR ">P1;frag_merged${count}\n";
    print PIR " : : : : : : : : : \n";
    print PIR "$pdb_seq_align*\n";
    close PIR;  

   `cp $outputfolder/frag_merge_dir/frag_merged$count.pir $outputfolder/frag_merge_dir/frag_merged$count/`;


  	open(RUNFILE,">$outputfolder/frag_merge_dir/frag_merged$count.sh") || die "Failed to write $outputfolder/frag_merge_dir/frag_merged$count.sh\n\n";
  	`touch $outputfolder/frag_merge_dir/frag_merged$count.sh.queued`;
  	print RUNFILE "#!/bin/bash\n\n";
  	print RUNFILE "mv $outputfolder/frag_merge_dir/frag_merged$count.sh.queued $outputfolder/frag_merge_dir/frag_merged$count.sh.running\n\n";
  	print RUNFILE "\nprintf \"Use MTMG to refine...\"\n";
   
    print RUNFILE "printf \"/storage/htc/bdm/jh7x3/Cryo_em_paper/paper_version_20190817/New_sequence_mapping_jiealgo/scripts/MTMG/mtmg $outputfolder/frag_merge_dir/frag_merged$count/ frag_merged$count.pir frag_merged$count  $outputfolder/frag_merge_dir/frag_merged$count/ /storage/htc/bdm/jh7x3/Cryo_em_paper/paper_version_20190817/New_sequence_mapping_jiealgo/scripts/MTMG/ /storage/htc/bdm/jh7x3/DeepRank/tools/R-3.2.0/bin/ 0 d\"\n\n";
  	print RUNFILE "/storage/htc/bdm/jh7x3/Cryo_em_paper/paper_version_20190817/New_sequence_mapping_jiealgo/scripts/MTMG/mtmg $outputfolder/frag_merge_dir/frag_merged$count/ frag_merged$count.pir frag_merged$count  $outputfolder/frag_merge_dir/frag_merged$count/ /storage/htc/bdm/jh7x3/Cryo_em_paper/paper_version_20190817/New_sequence_mapping_jiealgo/scripts/MTMG/ /storage/htc/bdm/jh7x3/DeepRank/tools/R-3.2.0/bin/ 0 d\n\n";
    print RUNFILE "cp $outputfolder/frag_merge_dir/frag_merged$count/frag_merged$count.pdb $outputfolder/frag_merge_dir/\n\n";
  	print RUNFILE "mv $outputfolder/frag_merge_dir/frag_merged$count.sh.running $outputfolder/frag_merge_dir/frag_merged$count.sh.done\n\n";
  	close RUNFILE;
  
    push @running_files,"$outputfolder/frag_merge_dir/frag_merged$count.sh";

}


MODELLING:
##### (7) Generate models

$shell_dir = "$outputfolder/frag_merge_dir/";

foreach $file_path (sort @running_files)
{
	## check the running jobs
	$min_elaps=0;
	while(1)
	{
		opendir(DIR,"$shell_dir") || die "Failed to open directory $shell_dir\n";
		@out_files = readdir(DIR);
		closedir(DIR);
		
		$running_num = 0;
		foreach $check_file (sort @out_files)
		{
			if($check_file eq '.' or $check_file eq '..' or substr($check_file,length($check_file)-8) ne '.running')
			{
				next;
			}
			$running_num++;
		}
		if($running_num<$proc_num)
		{
			last;
		}
		sleep(60);
		$min_elaps++;
		if($min_elaps > 60)
		{
			last; # move to next;
		}
	}
	
	if(!(-e substr($file_path,0,length($file_path)-3).".done"))
	{
		print "run test $file_path\n";
		system("sh $file_path &> $file_path.log &");
	}else{
		print "$file_path has been done\n";
		$queue_file = substr($file_path,0,length($file_path)-3).".queued";
		if(-e $queue_file)
		{
			`rm $queue_file`;
		}
	}
	
	$running_jobs=0;
	$processed_jobs=0;
	opendir(DIR,"$shell_dir") || die "Failed to open directory $shell_dir\n";
	@out_files = readdir(DIR);
	closedir(DIR);
	foreach $check_file (sort @out_files)
	{
		if($check_file eq '.' or $check_file eq '..')
		{
			next;
		}
		if(substr($check_file,length($check_file)-5) eq '.done')
		{
			$processed_jobs++;
		}
		if(substr($check_file,length($check_file)-8) eq '.running')
		{
			$running_jobs++;
		}
	}
	$remain_jobs = @running_files-$processed_jobs-$running_jobs;
	print "Current running jobs ($running_num), processed jobs ($processed_jobs), unprocessed jobs ($remain_jobs)\n\n";
	sleep(5);
}

#### check if all files have finished
print "#### check if all files have finished\n";

while(1)
{

	opendir(DIR,"$shell_dir") || die "Failed to open directory $shell_dir\n";
	@out_files = readdir(DIR);
	closedir(DIR);

  $running_num = 0;
  foreach $check_file (sort @out_files)
  {
  	if($check_file eq '.' or $check_file eq '..' or substr($check_file,length($check_file)-3) eq '.sh')
  	{
  		next;
  	}
   
    if(substr($check_file,length($check_file)-8) eq '.running' or substr($check_file,length($check_file)-7) eq '.queued')
    {
  	  $running_num++;
    }
  }
  
  if($running_num>0)
  {
    print "$running_num jobs are still running, please wait\n";
  }else{
    print "All training jobs are done\n\n";
    last;
  }
  
  sleep(60*5);
  
}



sub generate_gaps
{
	my $gnum = $_[0]; 	
	my $gaps = "";
	my $i;
	for ($i = 0; $i < $gnum; $i++)
	{
		$gaps .= "-"; 
	}
	return $gaps; 
}

sub generate_aa
{
	my $gnum = $_[0]; 	
	my $gaps = "";
	my $i;
	for ($i = 0; $i < $gnum; $i++)
	{
		$gaps .= "G"; 
	}
	return $gaps; 
}


sub parse_pdb_row{
	my $row = shift;
	my $param = shift;
	my $result;
	$result = substr($row,6,5) if ($param eq "anum");
	$result = substr($row,12,4) if ($param eq "aname");
	$result = substr($row,16,1) if ($param eq "altloc");
	$result = substr($row,17,3) if ($param eq "rname");
	$result = substr($row,22,5) if ($param eq "rnum");
	$result = substr($row,21,1) if ($param eq "chain");
	$result = substr($row,30,8) if ($param eq "x");
	$result = substr($row,38,8) if ($param eq "y");
	$result = substr($row,46,8) if ($param eq "z");
	die "Invalid row[$row] or parameter[$param]" if (not defined $result);
	$result =~ s/\s+//g;
	return $result;
}


# Fischer-Krause ordered permutation generator

sub permute (&@) {
        my $code = shift;
        my @idx = 0..$#_;
        while ( $code->(@_[@idx]) ) {
                my $p = $#idx;
                --$p while $idx[$p-1] > $idx[$p];
                my $q = $p or return;
                push @idx, reverse splice @idx, $p;
                ++$q while $idx[$p-1] > $idx[$q];
                @idx[$p-1,$q]=@idx[$q,$p-1];
        }
}

