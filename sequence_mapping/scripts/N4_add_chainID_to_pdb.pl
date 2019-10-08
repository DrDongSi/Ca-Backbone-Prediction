##########################################################################################################
#                     Function about adding chain ID to the model        				                 #
#																										 #
#										Renzhi Cao  													 #
#																										 #
#									    10/4/2012														 #
#																										 #
#																										 #
#									Revised at 10/4/2012	                         					 #
#																										 #
##########################################################################################################
#! /usr/bin/perl -w
=pod
You may freely copy and distribute this document so long as the copyright is left intact. You may freely copy and post unaltered versions of this document in HTML and Postscript formats on a web site or ftp site. Lastly, if you do something injurious or stupid
because of this document, I don't want to know about it. Unless it's amusing.
=cut
 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
  if (@ARGV != 2 && @ARGV!=3)
    { # @ARGV used in scalar context = number of args
      print("This program tries to add chain ID to the protein, in default, we will add chain A\n");
	  print("Input is the path for the model, and the chain ID if you want to specific use, then the output model\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME  address_input_pdb address_output_model chain_ID(optional, indefault is A) \n");
          print "Revised at 10/5, don't disturbe the pdb format, no chomp(), the pdb has some space there.\n";
      print("\n********** example******\n");
	  print("perl $PROGRAM_NAME \n");
      print("perl $PROGRAM_NAME \n");
	  exit(1) ;
    }
 my $starttime = localtime();
 print "\n The time started at : $starttime.\n";
 my($chainID)="A";     # in default is A
 my($input_pdb)=$ARGV[0];
 my($output_pdb)=$ARGV[1];
 if(@ARGV==3)
 {# have specific chain ID
    $chainID=$ARGV[2];
 }
 -e $input_pdb || die "Input pdb is not exists!\n";
 if(!-s $output_pdb)
 {
      open (File, "&gt;$output_pdb");
	  chmod (0777, $output_pdb); 
      close (File);
 }
 else
 {
	 print "$output_pdb already exists!\n";
 }

##########################################################################################################
#              Function about openning a directory and processing the files                				 #
#																										 #
#										Renzhi Cao  													 #
#																										 #
#									    12/27/2011														 #
#																										 #
##########################################################################################################
my($line,$IN,$OUT,@tem_split);
$OUT=new FileHandle ">$output_pdb";
$IN=new FileHandle "$input_pdb";
defined($IN) || die "Cannot open input pdb : $input_pdb\n";
while(defined($line=<$IN>))
{
	#chomp($line); should delete, because of the pdb format
	#$line=~s/\s+$//;  # remove the windows character
    @tem_split=split(/\s+/,$line);
	if($tem_split[0] ne "ATOM")
	{# output this line directly 
		print $OUT $line;
	}
	else
	{# put the chainID here
		for(substr($line,21,1))
		{
			$_=$chainID;
		}
		
		print $OUT $line;
	}

}
$IN->close();
$OUT->close();
 my $endtime = localtime();
 print  "\nThe time ended at : $endtime.\n";
