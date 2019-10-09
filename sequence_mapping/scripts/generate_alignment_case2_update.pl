#!/usr/bin/perl -w
our %AA3TO1 = qw(ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V);
our %AA1TO3 = reverse %AA3TO1;


$installation_dir = '/home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/';
if (@ARGV < 4)
{
	die "Error: need at least four parameters: <path of Ca trace> <path of fasta sequence> <output-directory> <number of cpus> <previous fitted fragments, optional>.\n";
}
#### 
$input_pdb = shift @ARGV;
$fasta_file = shift @ARGV;
$outputfolder = shift @ARGV;
$proc_num = shift @ARGV;
$fitted_fragments = shift @ARGV;

-d $outputfolder || `mkdir $outputfolder`;
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

## loading previous fitted fragments if exists
%previous_fitted_fragments = ();
if(-e $fitted_fragments)
{
  open(TMP,"$fitted_fragments") || die "Failed to open file $fitted_fragments\n";
  while(<TMP>)
  {
    $line = $_;
    chomp $line;
    @tmp = split(/\s/,$line);
    $pdb_path = $tmp[0];
    $start = $tmp[1];
    $end = $tmp[2];
    if(-e $pdb_path)
    {
      print "Loading $pdb_path with range $start-$end\n\n";
      $previous_fitted_fragments{$pdb_path} = "$start-$end";
    }else{
      next;
    }
  }
  close TMP;
}




#####  (1) initialize sequence in fragment for forward

`cp $input_pdb $outputfolder/input.pdb`;
`$installation_dir/tools/MTMG/tools/pulchra_306/pulchra -c -s $outputfolder/input.pdb`;

$init_pdb = "$outputfolder/input.pdb";
if(!(-e "$outputfolder/input.rebuilt.pdb"))
{
  $init_pdb = "$outputfolder/input.pdb";
}else{

  `$installation_dir/tools/MTMG/Scwrl4 -i $outputfolder/input.rebuilt.pdb -o  $outputfolder/input_scwrl.pdb`;
  if(!(-e "$outputfolder/input_scwrl.pdb"))
  {
    $init_pdb = "$outputfolder/input.rebuilt.pdb";
  }else{
    $init_pdb = "$outputfolder/input_scwrl.pdb";
  }
}

open OUTPDB, ">$outputfolder/temp0.pdb" or die "ERROR! Could not open $outputfolder/temp0.pdb\n";

open INPUTPDB, "$outputfolder/input_scwrl.pdb" or die "ERROR! Could not open $outputfolder/input_scwrl.pdb";
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

open OUTPDB, ">$outputfolder/input_r.pdb" or die "ERROR! Could not open $outputfolder/input_r.pdb\n";

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
  next if $atomtype ne 'CA';
  
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


`$installation_dir/tools/MTMG/tools/pulchra_306/pulchra -c -s $outputfolder/input_r.pdb`;

$init_pdb = "$outputfolder/input_r.pdb";
if(!(-e "$outputfolder/input_r.rebuilt.pdb"))
{
  $init_pdb = "$outputfolder/input_r.pdb";
}else{

  `$installation_dir/tools/MTMG/Scwrl4 -i $outputfolder/input_r.rebuilt.pdb -o  $outputfolder/input_r_scwrl.pdb`;
  if(!(-e "$outputfolder/input_r_scwrl.pdb"))
  {
    $init_pdb = "$outputfolder/input_r.rebuilt.pdb";
  }else{
    $init_pdb = "$outputfolder/input_r_scwrl.pdb";
  }
}

open OUTPDB, ">$outputfolder/temp0_r.pdb" or die "ERROR! Could not open $outputfolder/temp0_r.pdb\n";

open INPUTPDB, "$outputfolder/input_r_scwrl.pdb" or die "ERROR! Could not open $outputfolder/input_r_scwrl.pdb";
@lines_PDB = <INPUTPDB>;
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





##### (3) rebuild side-chain according to new aa

$init_pdb = "$outputfolder/temp0.pdb";

`$installation_dir/tools/MTMG/Scwrl4 -i $outputfolder/temp0.pdb -o  $outputfolder/temp0_scwrl.pdb`;
if(!(-e "$outputfolder/temp0_scwrl.pdb"))
{
  $init_pdb = "$outputfolder/temp0.pdb";
}else{
  $init_pdb = "$outputfolder/temp0_scwrl.pdb";
}


$init_pdb_reverse = "$outputfolder/temp0_r.pdb";

`$installation_dir/tools/MTMG/Scwrl4 -i $outputfolder/temp0_r.pdb -o  $outputfolder/temp0_r_scwrl.pdb`;
if(!(-e "$outputfolder/temp0_r_scwrl.pdb"))
{
  $init_pdb_reverse = "$outputfolder/temp0_r.pdb";
}else{
  $init_pdb_reverse = "$outputfolder/temp0_r_scwrl.pdb";
}


open INPUTPDB, "$init_pdb" or die "ERROR! Could not open $init_pdb\n";

@lines_PDB = <INPUTPDB>; 
close INPUTPDB; 
$pdb_seq="";
$pdb_record = ();
$this_rchain="";
foreach (@lines_PDB) {
	next if $_ !~ m/^ATOM/;
	next unless (parse_pdb_row($_,"aname") eq "CA");
	$this_rchain = parse_pdb_row($_,"chain");
	$rname = $AA3TO1{parse_pdb_row($_,"rname")};
	$rnum = parse_pdb_row($_,"rnum");
  if(exists($pdb_record{"$this_rchain-$rname-$rnum"}))
  {
    next;
  }else{
    $pdb_record{"$this_rchain-$rname-$rnum"} = 1;
  }
  
  $pdb_seq .= $rname;
}

print "PDB_seq: ".length($pdb_seq)."\n";
print "fas_seq: ".length($qseq)."\n";
print "$pdb_seq\n$qseq\n";

if(length($pdb_seq) > length($qseq))
{
  die "Fragment sequence is larger than fasta sequence. This program will be called if the length of fragment is shorter than fasta sequence. Checking $pdb_seq\n$qseq\n\n";
}


##### (4) generate alignment

$n_gaps = length($pdb_seq) - length($qseq);

-d "$outputfolder/Alignments" || `mkdir $outputfolder/Alignments/`;
-d "$outputfolder/Atoms" || `mkdir $outputfolder/Atoms/`;
-d "$outputfolder/Models" || `mkdir $outputfolder/Models/`;

`cp $init_pdb $outputfolder/Atoms/aaaaA.atm`;
`cp $init_pdb_reverse $outputfolder/Atoms/aaaaB.atm`;
chdir("$outputfolder/Atoms");
#`gzip -f aaaaA.atom`; 
#`gzip -f aaaaB.atom`; 


##### (6) Generate alignments for forward direction
chdir($outputfolder);


$shell_dir = "$outputfolder/sh_src";
-d $shell_dir || `mkdir $shell_dir`;

$n_gaps = length($qseq) - length($pdb_seq);
@running_files = ();


## insert gap in the N-terminal and C-terminal using forward direction

for($i=1;$i<=$n_gaps;$i++)
{
    $newseq = $pdb_seq;
    
    $start = $i;
    $end = $i + length($pdb_seq);
    
    ## check if the cover overlaps previous fitted structures
    $valid_pos = 1;
    foreach $fiited_pdb (sort keys %previous_fitted_fragments)
    {
      $fitted_info = $previous_fitted_fragments{$fiited_pdb};
      @tmp = split('-',$fitted_info);
      $fitted_start = $tmp[0];
      $fitted_end = $tmp[1];
      if(($start>=$fitted_start and  $start <= $fitted_end) or ($end>=$fitted_start and  $end <= $fitted_end))
      {
        $valid_pos = 0;
      }
      
    }
    if($valid_pos == 0)
    {
      next;
    }
    
    
    -d "$outputfolder/Alignments/temp_NC_$i/" || `mkdir $outputfolder/Alignments/temp_NC_$i/`;
    `cp $outputfolder/Atoms/aaaaA.atm $outputfolder/Alignments/temp_NC_$i/`;
    

    ### update the coordinates in previous fitting iterations 
      open INPUTPDB, "$outputfolder/Alignments/temp_NC_$i/aaaaA.atm" or die "ERROR! Could not open $outputfolder/Alignments/temp_NC_$i/aaaaA.atm\n";
      @lines_PDB = <INPUTPDB>; 
      close INPUTPDB; 
      %fitted_pdb_record = ();
      
      foreach (@lines_PDB) {
      	next if $_ !~ m/^ATOM/;
      	next unless (parse_pdb_row($_,"aname") eq "CA");
      	$this_rchain = parse_pdb_row($_,"chain");
      	$rname = $AA3TO1{parse_pdb_row($_,"rname")};
      	$atomCounter = parse_pdb_row($_,"anum");
      	$resCounter = $start + parse_pdb_row($_,"rnum")-1;

        $rnum_string = sprintf("%4s", $resCounter);
        $anum_string = sprintf("%5s", $atomCounter);
        $row = substr($_,0,6).$anum_string.substr($_,11,5)." ".substr($_,17,3)." "." ".$rnum_string." ".substr($_,27);        
        if(exists($fitted_pdb_record{"$resCounter"}))
        {
          next;
        }else{
          $fitted_pdb_record{"$resCounter"} = "$row";
          #print "$resCounter\n";
        }
      }

    $newseq = &generate_gaps($i).$newseq.&generate_gaps($n_gaps-$i);    
    #print "$newseq\n";
    
    
    
    ### loading previous fragments
    $fragidx = 0;
    foreach $fiited_pdb (sort keys %previous_fitted_fragments)
    {
      $fitted_info = $previous_fitted_fragments{$fiited_pdb};
      @tmp = split('-',$fitted_info);
      $fitted_start = $tmp[0];
      $fitted_end = $tmp[1];
      ## get atoms within ranges
      $fitted_pdb_seq = "";
       open INPUTPDB, "$fiited_pdb" or die "ERROR! Could not open $fiited_pdb\n";
      @lines_PDB = <INPUTPDB>; 
      close INPUTPDB; 
      
      foreach (@lines_PDB) {
      	next if $_ !~ m/^ATOM/;
      	next unless (parse_pdb_row($_,"aname") eq "CA");
      	$this_rchain = parse_pdb_row($_,"chain");
      	$this_x = parse_pdb_row($_,"x");
      	$this_y = parse_pdb_row($_,"y");
      	$this_z = parse_pdb_row($_,"z");
      	$rname = $AA3TO1{parse_pdb_row($_,"rname")};
      	$rnum = parse_pdb_row($_,"rnum");
        if($rnum>=$fitted_start and $rnum <= $fitted_end)
        {
          $fitted_pdb_seq .= $rname;
          if(exists($fitted_pdb_record{"$rnum"}))
          {
            next;
          }else{
            $fitted_pdb_record{"$rnum"} = $_;
          }
        }
      }
      
      $newseq = substr($newseq,0,$fitted_start).$fitted_pdb_seq.substr($newseq,$fitted_end+1); 
      #print "$newseq\n";   
      
    }
    #print "$qseq\n\n\n";
    
    open(OUTPUTPDB,">$outputfolder/Alignments/temp_NC_$i/ffffA.pdb") || die "Failed to write $outputfolder/Alignments/temp_NC_$i/ffffA.pdb\n";
    $min_rnum = 1000;
    $max_rnum = 0;
    foreach $rnum (sort {$a <=> $b} keys %fitted_pdb_record) {
      print OUTPUTPDB $fitted_pdb_record{$rnum};
      if($rnum <= $min_rnum )
      {
        $min_rnum = $rnum;
      }
      if($rnum >= $max_rnum )
      {
        $max_rnum = $rnum;
      }
    }
    close OUTPUTPDB;


    `$installation_dir/tools/MTMG/tools/pulchra_306/pulchra -c -s $outputfolder/Alignments/temp_NC_$i/ffffA.pdb`;
    
    $init_pdb = "$outputfolder/Alignments/temp_NC_$i/ffffA.pdb";
    if(!(-e "$outputfolder/Alignments/temp_NC_$i/ffffA.rebuilt.pdb"))
    {
      $init_pdb = "$outputfolder/Alignments/temp_NC_$i/ffffA.pdb";
    }else{
    
      `$installation_dir/tools/MTMG/Scwrl4 -i $outputfolder/Alignments/temp_NC_$i/ffffA.rebuilt.pdb -o  $outputfolder/Alignments/temp_NC_$i/ffffA_scwrl.pdb`;
      if(!(-e "$outputfolder/Alignments/temp_NC_$i/ffffA_scwrl.pdb"))
      {
        $init_pdb = "$outputfolder/Alignments/temp_NC_$i/ffffA.rebuilt.pdb";
      }else{
        $init_pdb = "$outputfolder/Alignments/temp_NC_$i/ffffA_scwrl.pdb";
      }
    }
    
    `cp $init_pdb $outputfolder/Alignments/temp_NC_$i/ffffA.atm`;



    
    $pir_file = "$outputfolder/Alignments/temp_NC_$i/align.pir";
    open(PIR, ">$pir_file") || die "can't create pir file $pir_file.\n";
    print PIR "C;cover size:X; local alignment length=X (original info = ffffA	X	X	0	0)\n";
    print PIR ">P1;ffffA\n";
    print PIR "structureX:ffffA: $min_rnum: $this_rchain: $max_rnum: : : : : \n"; 
    print PIR "$newseq*\n\n";
    
    
    print PIR "C;query_length:X num_of_temps:X cover_ratio:X cover:X not_cover:X\n"; 
    print PIR ">P1;temp_NC_$i\n";
    print PIR " : : : : : : : : : \n";
    print PIR "$qseq*\n";
    close PIR;  
    
  
  	open(RUNFILE,">$shell_dir/temp_NC_$i.sh") || die "Failed to write $shell_dir/temp_NC_$i.sh\n\n";
  	`touch $shell_dir/temp_NC_$i.sh.queued`;
  	print RUNFILE "#!/bin/bash\n\n";
  	print RUNFILE "mv $shell_dir/temp_NC_$i.sh.queued $shell_dir/temp_NC_$i.sh.running\n\n";
  	print RUNFILE "\nprintf \"Use MTMG to refine...\"\n";
  
  	print RUNFILE "mkdir -p $outputfolder/Alignments/temp_NC_$i/mtmg_refine/\n\n";
    print RUNFILE "printf \"$installation_dir/tools/MTMG/mtmg  $outputfolder/Alignments/temp_NC_$i/ align.pir temp_NC_$i $outputfolder/Alignments/temp_NC_$i/mtmg_refine/ $installation_dir/tools/MTMG/ $installation_dir/tools/R-3.2.0/bin/ 0 d\"\n\n";
  	print RUNFILE "$installation_dir/tools/MTMG/mtmg  $outputfolder/Alignments/temp_NC_$i/ align.pir temp_NC_$i $outputfolder/Alignments/temp_NC_$i/mtmg_refine/ $installation_dir/tools/MTMG/ $installation_dir/tools/R-3.2.0/bin/ 0 d\n\n";
    print RUNFILE "cp $outputfolder/Alignments/temp_NC_$i/mtmg_refine/temp_NC_$i.pdb $outputfolder/Models\n\n";
  
  	print RUNFILE "mv $shell_dir/temp_NC_$i.sh.running $shell_dir/temp_NC_$i.sh.done\n\n";
  	close RUNFILE;
  
    push @running_files,"$shell_dir/temp_NC_$i.sh";

}

## insert gap in the N-terminal and C-terminal using backward direction

for($i=1;$i<=$n_gaps;$i++)
{
    $newseq = $pdb_seq;
    
    $start = $i;
    $end = $i + length($pdb_seq);
    
    ## check if the cover overlaps previous fitted structures
    $valid_pos = 1;
    foreach $fiited_pdb (sort keys %previous_fitted_fragments)
    {
      $fitted_info = $previous_fitted_fragments{$fiited_pdb};
      @tmp = split('-',$fitted_info);
      $fitted_start = $tmp[0];
      $fitted_end = $tmp[1];
      if(($start>=$fitted_start and  $start <= $fitted_end) or ($end>=$fitted_start and  $end <= $fitted_end))
      {
        $valid_pos = 0;
      }
      
    }
    if($valid_pos == 0)
    {
      next;
    }
    
    $newseq = &generate_gaps($i).$newseq.&generate_gaps($n_gaps-$i);
    print "$newseq\n";
    
    -d "$outputfolder/Alignments/temp_r_NC_$i/" || `mkdir $outputfolder/Alignments/temp_r_NC_$i/`;
    `cp $outputfolder/Atoms/aaaaB.atm $outputfolder/Alignments/temp_r_NC_$i/`;
    
    

    ### update the coordinates in previous fitting iterations 
      open INPUTPDB, "$outputfolder/Alignments/temp_r_NC_$i/aaaaB.atm" or die "ERROR! Could not open $outputfolder/Alignments/temp_r_NC_$i/aaaaB.atm\n";
      @lines_PDB = <INPUTPDB>; 
      close INPUTPDB; 
      %fitted_pdb_record = ();
      
      foreach (@lines_PDB) {
      	next if $_ !~ m/^ATOM/;
      	next unless (parse_pdb_row($_,"aname") eq "CA");
      	$this_rchain = parse_pdb_row($_,"chain");
      	$rname = $AA3TO1{parse_pdb_row($_,"rname")};
      	$atomCounter = parse_pdb_row($_,"anum");
      	$resCounter = $start + parse_pdb_row($_,"rnum")-1;

        $rnum_string = sprintf("%4s", $resCounter);
        $anum_string = sprintf("%5s", $atomCounter);
        $row = substr($_,0,6).$anum_string.substr($_,11,5)." ".substr($_,17,3)." "." ".$rnum_string." ".substr($_,27);        
        if(exists($fitted_pdb_record{"$resCounter"}))
        {
          next;
        }else{
          $fitted_pdb_record{"$resCounter"} = "$row";
          #print "$resCounter\n";
        }
      } 
    
    
    
    ### loading previous fragments
    $fragidx = 0;
    foreach $fiited_pdb (sort keys %previous_fitted_fragments)
    {
      $fitted_info = $previous_fitted_fragments{$fiited_pdb};
      @tmp = split('-',$fitted_info);
      $fitted_start = $tmp[0];
      $fitted_end = $tmp[1];
      ## get atoms within ranges
      $fitted_pdb_seq = "";
       open INPUTPDB, "$fiited_pdb" or die "ERROR! Could not open $fiited_pdb\n";
      @lines_PDB = <INPUTPDB>; 
      close INPUTPDB; 
      
      foreach (@lines_PDB) {
      	next if $_ !~ m/^ATOM/;
      	next unless (parse_pdb_row($_,"aname") eq "CA");
      	$this_rchain = parse_pdb_row($_,"chain");
      	$this_x = parse_pdb_row($_,"x");
      	$this_y = parse_pdb_row($_,"y");
      	$this_z = parse_pdb_row($_,"z");
      	$rname = $AA3TO1{parse_pdb_row($_,"rname")};
      	$rnum = parse_pdb_row($_,"rnum");
        if($rnum>=$fitted_start and $rnum <= $fitted_end)
        {
          $fitted_pdb_seq .= $rname;
          if(exists($fitted_pdb_record{"$rnum"}))
          {
            next;
          }else{
            $fitted_pdb_record{"$rnum"} = $_;
          }
        }
      }
      
      $newseq = substr($newseq,0,$fitted_start).$fitted_pdb_seq.substr($newseq,$fitted_end+1); 
      print "$newseq\n";   
      
    }
    print "$qseq\n\n\n";
    
    open(OUTPUTPDB,">$outputfolder/Alignments/temp_r_NC_$i/ffffB.pdb") || die "Failed to write $outputfolder/Alignments/temp_r_NC_$i/ffffB.pdb\n";
    $min_rnum = 1000;
    $max_rnum = 0;
    foreach $rnum (sort {$a <=> $b} keys %fitted_pdb_record) {
      print OUTPUTPDB $fitted_pdb_record{$rnum};
      if($rnum <= $min_rnum )
      {
        $min_rnum = $rnum;
      }
      if($rnum >= $max_rnum )
      {
        $max_rnum = $rnum;
      }
    }
    close OUTPUTPDB;


    `$installation_dir/tools/MTMG/tools/pulchra_306/pulchra -c -s $outputfolder/Alignments/temp_r_NC_$i/ffffB.pdb`;
    
    $init_pdb = "$outputfolder/Alignments/temp_r_NC_$i/ffffB.pdb";
    if(!(-e "$outputfolder/Alignments/temp_r_NC_$i/ffffB.rebuilt.pdb"))
    {
      $init_pdb = "$outputfolder/Alignments/temp_r_NC_$i/ffffB.pdb";
    }else{
    
      `$installation_dir/tools/MTMG/Scwrl4 -i $outputfolder/Alignments/temp_r_NC_$i/ffffB.rebuilt.pdb -o  $outputfolder/Alignments/temp_r_NC_$i/ffffB_scwrl.pdb`;
      if(!(-e "$outputfolder/Alignments/temp_r_NC_$i/ffffB_scwrl.pdb"))
      {
        $init_pdb = "$outputfolder/Alignments/temp_r_NC_$i/ffffB.rebuilt.pdb";
      }else{
        $init_pdb = "$outputfolder/Alignments/temp_r_NC_$i/ffffB_scwrl.pdb";
      }
    }
    
    `cp $init_pdb $outputfolder/Alignments/temp_r_NC_$i/ffffB.atm`;



    
    $pir_file = "$outputfolder/Alignments/temp_r_NC_$i/align.pir";
    open(PIR, ">$pir_file") || die "can't create pir file $pir_file.\n";
    print PIR "C;cover size:X; local alignment length=X (original info = ffffB	X	X	0	0)\n";
    print PIR ">P1;ffffB\n";
    print PIR "structureX:ffffB: $min_rnum: $this_rchain: $max_rnum: : : : : \n"; 
    print PIR "$newseq*\n\n";
   
    print PIR "C;query_length:X num_of_temps:X cover_ratio:X cover:X not_cover:X\n"; 
    print PIR ">P1;temp_r_NC_$i\n";
    print PIR " : : : : : : : : : \n";
    print PIR "$qseq*\n";
    close PIR;  
    
  	open(RUNFILE,">$shell_dir/temp_r_NC_$i.sh") || die "Failed to write $shell_dir/temp_r_NC_$i.sh\n\n";
  	`touch $shell_dir/temp_r_NC_$i.sh.queued`;
  	print RUNFILE "#!/bin/bash\n\n";
  	print RUNFILE "mv $shell_dir/temp_r_NC_$i.sh.queued $shell_dir/temp_r_NC_$i.sh.running\n\n";
  	print RUNFILE "\nprintf \"Use MTMG to refine...\"\n";
  
  	print RUNFILE "mkdir -p $outputfolder/Alignments/temp_r_NC_$i/mtmg_refine/\n\n";
    print RUNFILE "printf \"$installation_dir/tools/MTMG/mtmg  $outputfolder/Alignments/temp_r_NC_$i/ align.pir temp_r_NC_$i $outputfolder/Alignments/temp_r_NC_$i/mtmg_refine/ $installation_dir/tools/MTMG/ $installation_dir/tools/R-3.2.0/bin/ 0 d\"\n\n";
  	print RUNFILE "$installation_dir/tools/MTMG/mtmg  $outputfolder/Alignments/temp_r_NC_$i/ align.pir temp_r_NC_$i $outputfolder/Alignments/temp_r_NC_$i/mtmg_refine/ $installation_dir/tools/MTMG/ $installation_dir/tools/R-3.2.0/bin/ 0 d\n\n";
    print RUNFILE "cp $outputfolder/Alignments/temp_r_NC_$i/mtmg_refine/temp_r_NC_$i.pdb $outputfolder/Models\n\n";
  	print RUNFILE "mv $shell_dir/temp_r_NC_$i.sh.running $shell_dir/temp_r_NC_$i.sh.done\n\n";
  	close RUNFILE;
  
    push @running_files,"$shell_dir/temp_r_NC_$i.sh";

}


##### (7) Generate models


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
	
	if(!(-e "$file_path.done"))
	{
		print "run test $file_path\n";
		system("sh $file_path &> $file_path.log &");
	}else{
		print "$file_path has been done\n";
		$queue_file = "$file_path.queued";
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







###################  convert coord back to fragment


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
