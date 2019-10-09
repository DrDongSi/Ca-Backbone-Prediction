#!/usr/bin/perl -w
 use FileHandle; # use FileHandles instead of open(),close()
 use Cwd;
 use Cwd 'abs_path';


######################## !!! Don't Change the code below##############

$install_dir = getcwd;
$install_dir=abs_path($install_dir);


if(!-s $install_dir)
{
	die "The CaTrace2Seq directory ($install_dir) is not existing, please revise the customize settings part inside the configure.pl, set the path as  your unzipped CaTrace2Seq directory\n";
}

if ( substr($install_dir, length($install_dir) - 1, 1) ne "/" )
{
        $install_dir .= "/";
}


print "checking whether the configuration file run in the installation folder ...";
$cur_dir = `pwd`;
chomp $cur_dir;
$configure_file = "$cur_dir/setup_env.pl";
if (! -f $configure_file || $install_dir ne "$cur_dir/")
{
        die "\nPlease check the installation directory setting and run the configure program under the main directory of CaTrace2Seq.\n";
}
print " OK!\n";


$CaTrace2Seq_db_tools_dir = "$install_dir/tools/";

if(!(-d $CaTrace2Seq_db_tools_dir))
{
	$status = system("mkdir $CaTrace2Seq_db_tools_dir");
	if($status)
	{
		die "Failed to create folder $CaTrace2Seq_db_tools_dir\n\n";
	}
}
$CaTrace2Seq_db_tools_dir=abs_path($CaTrace2Seq_db_tools_dir);



if ( substr($CaTrace2Seq_db_tools_dir, length($CaTrace2Seq_db_tools_dir) - 1, 1) ne "/" )
{
        $CaTrace2Seq_db_tools_dir .= "/";
}

print "Start install CaTrace2Seq into <$CaTrace2Seq_db_tools_dir>\n"; 



chdir($CaTrace2Seq_db_tools_dir);

$tools_dir = "$CaTrace2Seq_db_tools_dir";


if(!-d $tools_dir)
{ 
	$status = system("mkdir $tools_dir");
	if($status)
	{
		die "Failed to create folder ($tools_dir), check permission or folder path\n";
	}
	`chmod -R 755 $tools_dir`;
}




print "#########  (1) Configuring option files\n";

$option_list = "$install_dir/installation/configure_list";

if (! -f $option_list)
{
        die "\nOption file $option_list not exists.\n";
}
configure_file2($option_list,'bin');
configure_file2($option_list,'scripts');
configure_file2($option_list,'example');

print "#########  Configuring option files, done\n\n\n";

system("chmod +x $install_dir/bin/*.sh");
system("chmod +x $install_dir/examples/*.sh");




#### (2) Download basic tools
print("\n#### (2) Download basic tools\n\n");

chdir($tools_dir);
$basic_tools_list = "TMscore.tar.gz;qprob_package.tar.gz;R-3.2.0.tar.gz";
#$basic_tools_list = "TMscore.tar.gz;R-3.2.0.tar.gz";
@basic_tools = split(';',$basic_tools_list);
foreach $tool (@basic_tools)
{
	$toolname = substr($tool,0,index($tool,'.tar.gz'));
	if(-d "$tools_dir/$toolname")
	{
		if(-e "$tools_dir/$toolname/download.done")
		{
			print "\t$toolname is done!\n";
			next;
		}
	}elsif(-f "$tools_dir/$toolname")
	{
			print "\t$toolname is done!\n";
			next;
	}
	if(-e $tool)
	{
		 `rm $tool`;
	}
	if($tool eq 'R-3.2.0.tar.gz')
	{
		`wget http://sysbio.rnet.missouri.edu/bdm_download/DeepRank_db_tools/tools/$tool`;
	}else{
		`wget http://sysbio.rnet.missouri.edu/multicom_db_tools/tools/$tool`;
	}
	
	if(-e "$tool")
	{
		print "\n\t$tool is found, start extracting files......\n\n";
		`tar -zxf $tool`;
		if(-d $toolname)
		{
			`echo 'done' > $toolname/download.done`;
		}
		`rm $tool`;
		`chmod -R 755 $toolname`;
	}else{
		die "Failed to download $tool from http://sysbio.rnet.missouri.edu/multicom_db_tools/tools, please contact chengji\@missouri.edu\n";
	}
}

$tooldir = $CaTrace2Seq_db_tools_dir.'/qprob_package';
if(-d $tooldir)
{
	print "\n\n#########  Setting up qprob_package/\n";
	chdir $tooldir;
	if(-f 'configure.pl')
	{
		$status = system("perl configure.pl 2>&1 &> /dev/null");
		if($status){
			die "Failed to run perl configure.pl, possible reason is the permission conflict or incorrect software installation.\nIf the database and tools have already been configured, repeated configuration is not necessary.\n";
			exit(-1);
		}
	}else{
		die "The configure.pl file for $tooldir doesn't exist, please contact us(Jie Hou: jh7x3\@mail.missouri.edu)\n";
	}
	
	#### check the pspro and sspro in qprob
	print "#### check the pspro and sspro in qprob\n\n";
	for($num = 1; $num <=2; $num++)
	{
		$check_pass = 1;
		`$install_dir/tools/qprob_package/tools/pspro2.0/server/predict_seq_ss &> $install_dir/tools/qprob_package/tools/pspro2.0/server/test.log`;

		open(TMP,"$install_dir/tools/qprob_package/tools/pspro2.0/server/test.log");
		@info_log = <TMP>;
		close TMP;

		$infos = shift @info_log;

		if(index($infos,'model_definition')>0 or index($infos,'dataset_file') >0)
		{
			print "pspro is working for qprob\n";	
		}else{
			print "pspro is not working, checking backup files\n";
			`$install_dir/tools/pspro2_server_32bit/predict_seq_ss &> $install_dir/tools/pspro2_server_32bit/test.log`;

			open(TMP,"$install_dir/tools/pspro2_server_32bit/test.log");
			@info_log2 = <TMP>;
			close TMP;

			$infos2 = shift @info_log2;	
			if(index($infos2,'model_definition')>0 or index($infos2,'dataset_file') >0)
			{
				print "backup pspro is working\n";	
				`cp -ar $install_dir/tools/pspro2_server_32bit/* $install_dir/tools/qprob_package/tools/pspro2.0/server/`;
			}else{	
				print "Both pspro binary versions failed to pass the examination, qprob may not work\n";
			}
			$check_pass = 0;
		}


	
		`$install_dir/tools/qprob_package/tools/sspro4/server/predict_seq_ss &> $install_dir/tools/qprob_package/tools/sspro4/server/test.log`;

		open(TMP,"$install_dir/tools/qprob_package/tools/sspro4/server/test.log");
		@info_log = <TMP>;
		close TMP;

		$infos = shift @info_log;

		if(index($infos,'model_definition')>0 or index($infos,'dataset_file') >0)
		{
			print "sspro is working for qprob\n";	
		}else{
			print "sspro is not working, checking backup files\n";
			`$install_dir/tools/pspro2_server_32bit/predict_seq_ss &> $install_dir/tools/pspro2_server_32bit/test.log`;

			open(TMP,"$install_dir/tools/pspro2_server_32bit/test.log");
			@info_log2 = <TMP>;
			close TMP;

			$infos2 = shift @info_log2;	
			if(index($infos2,'model_definition')>0 or index($infos2,'dataset_file') >0)
			{
				print "backup pspro is working\n";	
				`cp -ar $install_dir/tools/pspro2_server_32bit/* $install_dir/tools/qprob_package/tools/sspro4/server/`;
			}else{	
				print "Both pspro binary versions failed to pass the examination, qprob may not work\n";
			}
			$check_pass = 0;
		}
		if($check_pass ==1)
		{
			last; # pass 
		}

	}
}

$addr_scwrl4 = $CaTrace2Seq_db_tools_dir."/MTMG/";
if(-d $addr_scwrl4)
{
	print "\n#########  Setting up scwrl4 \n";
	$addr_scwrl_orig = $addr_scwrl4."/"."Scwrl4.ini";
	$addr_scwrl_back = $addr_scwrl4."/"."Scwrl4.ini.back";
	system("cp $addr_scwrl_orig $addr_scwrl_back");
	@ttt = ();
	$OUT = new FileHandle ">$addr_scwrl_orig";
	$IN=new FileHandle "$addr_scwrl_back";
	while(defined($line=<$IN>))
	{
		chomp($line);
		@ttt = split(/\s+/,$line);
		
		if(@ttt>1 && $ttt[1] eq "FilePath")
		{
			print $OUT "\tFilePath\t=\t$addr_scwrl4/bbDepRotLib.bin\n"; 
		}
		else
		{
			print $OUT $line."\n";
		}
	}
	$IN->close();
	$OUT->close();
	print "Done\n";
}


print "\n\n";


#### install R-3.2.0.tar.gz

open(OUT,">$install_dir/installation/P1_install_R-3.2.0.sh") || die "Failed to open file $install_dir/installation/P1_install_R-3.2.0.sh";
print OUT "#!/bin/bash -e\n\n";
print OUT "echo \" Start compile R-3.2.0 (will take ~3 min)\"\n\n";
print OUT "cd $install_dir/tools/R-3.2.0\n\n";
print OUT "make clean\n\n";
print OUT "./configure --prefix=$install_dir/tools/R-3.2.0  --with-readline=no --with-x=no\n\n";
print OUT "make\n\n";
print OUT "make install\n\n";
print OUT "echo \"installed\" > $install_dir/tools/R-3.2.0/install.done\n\n";
close OUT;




sub prompt_yn {
  my ($query) = @_;
  my $answer = prompt("$query (Y/N): ");
  return lc($answer) eq 'y';
}
sub prompt {
  my ($query) = @_; # take a prompt string as argument
  local $| = 1; # activate autoflush to immediately show the prompt
  print $query;
  chomp(my $answer = <STDIN>);
  return $answer;
}


sub configure_file{
	my ($option_list,$prefix) = @_;
	open(IN,$option_list) || die "Failed to open file $option_list\n";
	$file_indx=0;
	while(<IN>)
	{
		$file = $_;
		chomp $file;
		if ($file =~ /^$prefix/)
		{
			$option_default = $install_dir.$file.'.default';
			$option_new = $install_dir.$file;
			$file_indx++;
			print "$file_indx: Configuring $option_new\n";
			if (! -f $option_default)
			{
					die "\nOption file $option_default not exists.\n";
			}	
			
			open(IN1,$option_default) || die "Failed to open file $option_default\n";
			open(OUT1,">$option_new") || die "Failed to open file $option_new\n";
			while(<IN1>)
			{
				$line = $_;
				chomp $line;

				if(index($line,'SOFTWARE_PATH')>=0)
				{
					$line =~ s/SOFTWARE_PATH/$install_dir/g;
					$line =~ s/\/\//\//g;
					print OUT1 $line."\n";
				}else{
					print OUT1 $line."\n";
				}
			}
			close IN1;
			close OUT1;
		}
	}
	close IN;
}


sub configure_tools{
	my ($option_list,$prefix,$DBtool_path) = @_;
	open(IN,$option_list) || die "Failed to open file $option_list\n";
	$file_indx=0;
	while(<IN>)
	{
		$file = $_;
		chomp $file;
		if ($file =~ /^$prefix/)
		{
			$option_default = $DBtool_path.$file.'.default';
			$option_new = $DBtool_path.$file;
			$file_indx++;
			print "$file_indx: Configuring $option_new\n";
			if (! -f $option_default)
			{
					next;
					#die "\nOption file $option_default not exists.\n";
			}	
			
			open(IN1,$option_default) || die "Failed to open file $option_default\n";
			open(OUT1,">$option_new") || die "Failed to open file $option_new\n";
			while(<IN1>)
			{
				$line = $_;
				chomp $line;

				if(index($line,'SOFTWARE_PATH')>=0)
				{
					$line =~ s/SOFTWARE_PATH/$DBtool_path/g;
					$line =~ s/\/\//\//g;
					print OUT1 $line."\n";
				}else{
					print OUT1 $line."\n";
				}
			}
			close IN1;
			close OUT1;
		}
	}
	close IN;
}



sub configure_file2{
	my ($option_list,$prefix) = @_;
	open(IN,$option_list) || die "Failed to open file $option_list\n";
	$file_indx=0;
	while(<IN>)
	{
		$file = $_;
		chomp $file;
		if ($file =~ /^$prefix/)
		{
			@tmparr = split('/',$file);
			$filename = pop @tmparr;
			chomp $filename;
			$filepath = join('/',@tmparr);
			$option_default = $install_dir.$filepath.'/.'.$filename.'.default';
			$option_new = $install_dir.$file;
			$file_indx++;
			print "$file_indx: Configuring $option_new\n";
			if (! -f $option_default)
			{
					die "\nOption file $option_default not exists.\n";
			}	
			
			open(IN1,$option_default) || die "Failed to open file $option_default\n";
			open(OUT1,">$option_new") || die "Failed to open file $option_new\n";
			while(<IN1>)
			{
				$line = $_;
				chomp $line;

				if(index($line,'SOFTWARE_PATH')>=0)
				{
					$line =~ s/SOFTWARE_PATH/$install_dir/g;
					$line =~ s/\/\//\//g;
					print OUT1 $line."\n";
				}else{
					print OUT1 $line."\n";
				}
			}
			close IN1;
			close OUT1;
		}
	}
	close IN;
}



