BOOSTLIB = -lboost_regex-mt -lpthread

MYSQLLIB = -lmysqlclient -I/usr/include/mysql/

MYSQLPPLIB = -lmysqlpp -I/usr/include/mysql++/

CXX = icpc
CXXFLAGS = -w0 -fast
#CXX = g++
#CXXFLAGS = -O3

ifdef OPENMP
CXXFLAGS += -DOPENMP=1 -openmp
endif

CPP = $(CXX) $(CXXFLAGS)
RANLIB = ranlib
OBJEXT = o
AR = ar
ARFLAGS = cru
WNLIB_AR = $(AR) $(ARFLAGS)
WNLIB_OBJECTS = wn1dmin.$(OBJEXT) wnasrt.$(OBJEXT) \
	wnbar.$(OBJEXT) wnbar2.$(OBJEXT) \
	wncnjfg.$(OBJEXT) wnconjd.$(OBJEXT) wnconjg.$(OBJEXT) \
	wndcmp.$(OBJEXT) wndot.$(OBJEXT) wnffln.$(OBJEXT) \
	wnlbfgs.$(OBJEXT) wnmbtr.$(OBJEXT) wnmecp.$(OBJEXT) \
	wnmem.$(OBJEXT) wnmemb.$(OBJEXT) wnmemg.$(OBJEXT) \
	wnmemn.$(OBJEXT) wnmmk.$(OBJEXT) wnmqf.$(OBJEXT) \
	wnnlp.$(OBJEXT) wnnmgd.$(OBJEXT) wnparvect.$(OBJEXT) \
	wnpen.$(OBJEXT) wnpoly.$(OBJEXT) wnqfit.$(OBJEXT) \
	wnrdb.$(OBJEXT) wnrflt.$(OBJEXT) wnrib.$(OBJEXT) \
	wnrmod.$(OBJEXT) wnrnd.$(OBJEXT) wnrnrm.$(OBJEXT) \
	wnrnzsgn.$(OBJEXT) wnroot.$(OBJEXT) wnrtab.$(OBJEXT) \
	wnscnt.$(OBJEXT) wntol.$(OBJEXT) wntwp.$(OBJEXT) \
	wnvadd2.$(OBJEXT) wnvadd3.$(OBJEXT) wnvcpy.$(OBJEXT) \
	wnvgen.$(OBJEXT) wnvmk.$(OBJEXT) wnvmlt.$(OBJEXT) \
	wnvnrm.$(OBJEXT) wnvprn.$(OBJEXT) wnvsub2.$(OBJEXT)
WNLIB = libwnlib.a
WNLIB_SOURCES = wn1dmin.c wnasrt.c wnbar.c wnbar2.c\
	 wncnjfg.c wnconjd.c wnconjg.c\
	 wndcmp.c wndot.c wnffln.c\
	wnlbfgs.c wnmbtr.c wnmecp.c wnmem.c\
	 wnmemb.c wnmemg.c wnmemn.c  wnmmk.c wnmqf.c\
	 wnnlp.c wnnmgd.c wnparvect.c wnpen.c\
	wnpoly.c wnqfit.c wnrdb.c\
	wnrflt.c wnrib.c wnrmod.c wnrnd.c wnrnrm.c\
	wnrnzsgn.c wnroot.c wnrtab.c\
	wnscnt.c wntol.c wntwp.c wnvadd2.c\
	wnvadd3.c wnvcpy.c wnvgen.c wnvmk.c wnvmlt.c wnvnrm.c wnvprn.c\
	wnvsub2.c

kDP_SOURCES = kDP.cpp matrix.cpp sequence.cpp fasta_db_reader.cpp packed_sequence.cpp recycle_table.cpp kmer_seed_list.cpp kmer.cpp kmeraln.cpp pair_aln.cpp profile.cpp converter.cpp kmer_matrix.cpp
kDP_FLAGS = -D CM_NO_BOOST

kClust_SOURCES = kClust.cpp k_clust.cpp chronometer.cpp matrix.cpp sequence.cpp fasta_db_reader.cpp packed_sequence.cpp recycle_table.cpp base_table.cpp kmer_seed_list.cpp kmer.cpp kmeraln.cpp pair_aln.cpp params.cpp clusters.cpp block_manager.cpp profile.cpp converter.cpp kmer_matrix.cpp

kClust_mkAln_SOURCES = kClust_mkAln.cpp matrix.cpp sequence.cpp fasta_db_reader.cpp packed_sequence.cpp alignment.cpp mkaln.cpp clusters.cpp params.cpp profile.cpp converter.cpp  kmer_matrix.cpp

#kClust_mkAln_onlycons_SOURCES = kClust_mkAln.cpp matrix.cpp sequence.cpp fasta_db_reader.cpp packed_sequence.cpp alignment.cpp mkaln_onlycons.cpp clusters.cpp params.cpp profile.cpp converter.cpp  kmer_matrix.cpp

kClust_benchmark_SOURCES = kClust_benchmarks.cpp k_clust.cpp chronometer.cpp matrix.cpp sequence.cpp fasta_db_reader.cpp packed_sequence.cpp recycle_table.cpp base_table.cpp kmer_seed_list.cpp kmer.cpp kmeraln.cpp pair_aln.cpp params.cpp clusters.cpp block_manager.cpp profile.cpp converter.cpp kmer_matrix.cpp
 
TARGETS = kDP kClust kClust_mkAln kClust_benchmark

all: $(TARGETS)

libwnlib.a: $(WNLIB_OBJECTS) 
	-rm -f libwnlib.a
	$(WNLIB_AR) libwnlib.a $(WNLIB_OBJECTS)
	$(RANLIB) libwnlib.a

kDP: 	$(WNLIB) $(kDP_SOURCES)
	$(CPP) $(kDP_FLAGS) -o kDP $(kDP_SOURCES) $(WNLIB)


kClust: $(WNLIB) $(kClust_SOURCES) 
	$(CPP) -o kClust $(kClust_SOURCES) $(WNLIB) $(BOOSTLIB)

#kClust_mkAlnDB: $(WNLIB)
#	$(CPP) -o kClust_mkAlnDB $(kClust_mkAlnDB_SOURCES) $(WNLIB) $(BOOSTLIB) $(MYSQLLIB) $(MYSQLPPLIB)

kClust_mkAln: $(WNLIB) $(kClust_mkAln_SOURCES)
	$(CPP) -o kClust_mkAln $(kClust_mkAln_SOURCES) $(WNLIB) $(BOOSTLIB)

#kClust_mkAln_onlycons: $(WNLIB)
#	$(CPP) -o kClust_mkAln_onlycons $(kClust_mkAln_onlycons_SOURCES) $(WNLIB) $(BOOSTLIB)
	
kClust_benchmark: $(WNLIB) $(kClust_benchmark_SOURCES) 
	$(CPP) -o kClust_benchmark $(kClust_benchmark_SOURCES) $(WNLIB) $(BOOSTLIB)

.PHONY: all clean-wnlib clean

clean-wnlib: 
	-rm -f $(libwnlib_a_OBJECTS)
	-rm -f libwnlib.a

clean: clean-wnlib
	-rm -f *.o
	-rm -f $(TARGETS)
