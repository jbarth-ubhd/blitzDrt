#!/usr/bin/perl
printf "#define SOURCECODE_TIMESTAMP \"";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
                       localtime(time);
my @abbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
printf "%s @ %s %2d %4d %2d:%02d:%02d\"\n",
"cpio, gzipped", $abbr[$mon], $mday, $year+1900, $hour, $min, $sec;

printf "const unsigned char SOURCECODE[]=";
open my $f, "-|", "find . -depth -mount |grep -v productiveDrt |grep -v sourcecode.h |cpio -o --quiet |gzip";
my $pre="{";
my $bytes=0;
while(<$f>) {
	for my $char (split //, $_) {
		print $pre;
		if($bytes++%20==19) { print "\n"; }
		print ord($char);
		$pre=",";
	}
}
printf "};\n";
printf "#define SOURCECODE_LENGTH %d\n", $bytes;
