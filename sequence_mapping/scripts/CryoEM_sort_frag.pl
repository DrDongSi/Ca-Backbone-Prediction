##############  This program will be called if the length of fragment is longer than fasta sequence


use List::Util qw( min max );
our %AA3TO1 = qw(ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V);
our %AA1TO3 = reverse %AA3TO1;

sub permute (&@);  # predeclare sub, paste sub at bottom
 
if (@ARGV != 2)
{
	die "Error: need seven parameters: domain_list, domain model folder, query file(fasta), target id, output dir, modeller path, model number. \n";
}
#### 
$input_pdb = shift @ARGV;
$outputfolder = shift @ARGV;


-d $outputfolder || `mkdir $outputfolder`;

$fragseq = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

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
  #if($frag_len<10)
  #{
  #  next;
  #}
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

