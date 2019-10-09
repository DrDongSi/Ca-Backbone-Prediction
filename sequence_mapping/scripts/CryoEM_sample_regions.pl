##############  This program will be called if the length of fragment is longer than fasta sequence


use List::Util qw( min max );
our %AA3TO1 = qw(ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V);
our %AA1TO3 = reverse %AA3TO1;

sub permute (&@);  # predeclare sub, paste sub at bottom
 
if (@ARGV != 1)
{
	die "Error: need seven parameters: domain_list, domain model folder, query file(fasta), target id, output dir, modeller path, model number. \n";
}
#### 
$input_pdb = shift @ARGV;  #frag_merged1/aaaaA.atm

open INPUTPDB, $input_pdb or die "ERROR! Could not open $input_pdb";
@lines_PDB = <INPUTPDB>;
close INPUTPDB;


$prevrNum = "XX";
$num= 0;
%pdb_lines = ();
foreach (@lines_PDB) {
	$line = $_;
	chomp $line;
 
	next if $line !~ m/^ATOM/;
	@tmp = split(/\s+/,$line);
  $atomtype = parse_pdb_row($line,"aname");
  next if $atomtype != 'CA';
  $this_rnum = parse_pdb_row($line,"rnum");
  $pdb_lines{$this_rnum} = $line;
  $num++;
  if($num==1)
  {
    $prevrNum = $this_rnum;
    next;
  }
  if ($this_rnum - $prevrNum > 1 ) {
         ## get distance 
        $line1_first = $pdb_lines{$prevrNum};
        $line1_first_x = parse_pdb_row($line1_first,"x");
        $line1_first_y = parse_pdb_row($line1_first,"y");
        $line1_first_z = parse_pdb_row($line1_first,"z");
        
        $line2_first = $pdb_lines{$this_rnum};
        $line2_first_x = parse_pdb_row($line2_first,"x");
        $line2_first_y = parse_pdb_row($line2_first,"y");
        $line2_first_z = parse_pdb_row($line2_first,"z");
        
        $distance1 = sqrt(($line1_first_x-$line2_first_x)*($line1_first_x-$line2_first_x) + ($line1_first_y-$line2_first_y)*($line1_first_y-$line2_first_y)+($line1_first_z-$line2_first_z)*($line1_first_z-$line2_first_z));
        
        print "$prevrNum-$this_rnum:$distance1\n";
  }
  $prevrNum = $this_rnum;
  
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

