ACLOCAL_AMFLAGS = -I m4
SUBDIRS = libs igor_src

# modelsdir = $(pkgdatadir) 
# models_DATA = models demo

#Since automake cannot copy directly directories another workaround:

demodir = $(pkgdatadir)/demo
demo_DATA = demo/genomicVs_with_primers.fasta demo/genomicDs.fasta demo/genomicJs_all_curated.fasta demo/murugan_naive1_noncoding_demo_seqs.txt

install-data-local:
	cp -r models/ $(pkgdatadir)

uninstall-local:
	rm -rf $(pkgdatadir)/models 
