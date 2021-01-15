if [ $# -eq 0 ]; then
    echo "Syntax"
    echo "     $(basename $0)  File1.csv  File2.csv  File3.csv  ..."
    echo
    echo "Each of the files is a csv output file from Fluidigm in the format specified at"
    echo "http://161.55.237.67/~eriq/dokuwiki/doku.php?id=assist:clemento:conversion_input_format"
    echo "which is a web page accessible from within the lab."
    echo
    echo "The output is a series of files named:"
    echo "     File1_twocol.txt File2_twocol.txt ...      and"
    echo "     File1_gsisim.txt File2_gsisim.txt ..."
    echo "and it also creates two compiled files named"
    echo "     AllFiles_TwoCol.txt and AllFiles_gsisim.txt"
    echo "which contain all the individuals from all the files"
    echo "merged together."
    echo
    echo "This file is under version control.  On Eric's hard drive it exists at:"
    echo "/Users/eriq/Documents/xp_dev_svn_checkouts/gsi_sim/snpset/2010_SNPset_GSI_TOOLS/DataConversion/script/fdoCSV2tcf_and_gs.sh";
    Rev="SVN Revision Number"
    echo $Rev:$;

    exit 1;
fi


###################################################
# wrap up the awk script as a function that takes a filename as an argument
function MakeTwoCol
{
    # Intitialize the zero report file to empty
    rm -f ZeroReport.txt;
    touch ZeroReport.txt;

    cat $1 | awk '
    BEGIN {
    # make tab the field separator
	FS=",";

    # now define the code maps
    map["A:A"] = "1\t1";
    map["A:C"] = "1\t2";
    map["A:G"] = "1\t3";
    map["A:T"] = "1\t4";
    map["A:-"] = "1\t5";
    map["C:A"] = "2\t1";
    map["C:C"] = "2\t2";
    map["C:G"] = "2\t3";
    map["C:T"] = "2\t4";
    map["C:-"] = "2\t5";
    map["G:A"] = "3\t1";
    map["G:C"] = "3\t2";
    map["G:G"] = "3\t3";
    map["G:T"] = "3\t4";
    map["G:-"] = "3\t5";
    map["T:A"] = "4\t1";
    map["T:C"] = "4\t2";
    map["T:G"] = "4\t3";
    map["T:T"] = "4\t4";
    map["T:-"] = "4\t5";
    map["-:A"] = "5\t1";
    map["-:C"] = "5\t2";
    map["-:G"] = "5\t3";
    map["-:T"] = "5\t4";
    map["-:-"] = "5\t5";

    map["X:X"] = "6\t6";  # these are for the steelhead sex locus
    map["X:Y"] = "6\t7";
    map["Y:X"] = "7\t6";
    map["Y:Y"] = "7\t7";  # we should never see this.  I could omit it and have it produce a zero report finding.
}


# look for the tag that tells us to start processing lines
$1~/SNP Converted Calls/ {go=1; n=0; next}

# skip the line immediately after that
go==1 && n==0 {++n; next}

# duplicate the header names on the next line.  Note, locus names start at column 3
go==1 && n==1 {
    ++n;
    printf("Sample_ID");
    for(i=3;i<=NF;i++) printf("\t%s\t%s",$i,$i);
    printf("\n");
    next;
}


# now, on all the other lines we just rip through them, so long as they are not NTC (negative controls)
go==1 && n==2 && $2!="NTC" {
    printf("%s",$2);
    for(i=3;i<=NF;i++) {
	if( $i in map ) {
	    printf("\t%s",map[$i]);
	}
	else {
	    # here we keep track of what symbols we dole out zeroes for
	    zeros[$i]++;

	    # and here we print out the zeroes
	    printf("\t0\t0");
	}

    }
    printf("\n");
}




# then, at the end we produce our report of what was called a 0
END {
    for(i in zeros) print i,zeros[i] > "ZeroReport.txt"
}
'
}
# That is the end of that function MakeTwoCol
#############################################################




######################################################
# here is a function for converting two column to gsi_sim format
# you call it with a single string which is the name you want to
# give to the fishery.
function TwoCol2GsiSim {

    rm -f xxxtempxxx;
    awk -v fn=$1 '
     NR==1 {
     for(i=2;i<=NF;i+=2) print $i;
     print "POP",fn;
     next;
     }
     { print }
   ' $1  > xxxtempxxx

    NL=$( awk 'NF==1 {n++} END {print n}' xxxtempxxx);
    NF=$( awk 'NF>10 {n++} END {print n}' xxxtempxxx);

    (
	echo "$NF $NL";
	cat xxxtempxxx;
	)

}



# that is the end of TwoCol2GsiSim
###########################################################


# remove the AllFiles_TwoCol.txt file
rm -f AllFiles_TwoCol.txt;

times=0;

for i in "$@"; do

  # make new names for things
  j=$(basename $i)
  k=${j/.csv/_twocol.txt}
  l=${j/.csv/_gsisim.txt}

  # make the two col format for individual files
  MakeTwoCol $i > $k;
  mv ZeroReport.txt ZeroReport_${k}
  echo "Created $k from $i"

  # now, make the gsi_sim file
  TwoCol2GsiSim $k > $l;
  echo "   Created $l from $k"


  # if first file, copy it all to the AllFiles_TwoCol.txt file
  if [ $times -eq 0 ]; then
      cat $k > AllFiles_TwoCol.txt;
  else
      awk 'NR>1' $k >> AllFiles_TwoCol.txt;
  fi
  echo "   Added $k to AllFiles_TwoCol.txt"


  # increment the counter by 1
  times=$((times + 1))

done


# and, in the very end, we are going to create the AllFiles gsi_sim file:
TwoCol2GsiSim AllFiles_TwoCol.txt > AllFiles_gsisim.txt;

echo "Created  AllFiles_gsisim.txt from file AllFiles_TwoCol.txt"


# and at the end we glom the zero reports into a single file:
for i in ZeroReport*; do echo; echo $i; cat $i; done > AllZeroReports.txt;

echo "Glommed the ZeroReports into the file AllZeroReports.txt"

