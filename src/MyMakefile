# Jesse Windle, 2019

# R.home("include")
RINC := $(shell Rscript --slave -e "cat(R.home('include'))")
# R.home("lib")
RLIB := $(shell Rscript --slave -e "cat(R.home('lib'))")

CC = g++-8

OBJ = polyagamma_hybrid.o PolyaGammaApproxSP.o InvertY.o PolyaGammaApproxAlt.o PolyaGamma.o truncated_norm.o truncated_gamma.o inverse_gaussian.o

polyagamma.so: $(OBJ)
	$(CC) $(OBS) -I $(RINC) -L $(RLIB) polyagamma.so

%.o: %.cpp
	$(CC) $< -I $(RINC) -c -o $@


# PolyaGammaApproxSP.o: PolyaGammaApproxSP.h PolyaGammaApproxSP.cpp
# 	$(CC) PolyaGammaApproxSP.cpp -I $(RINC) -c -o PolyaGammaApproxSP.o

# InvertY.o: InvertY.h InvertY.cpp
# 	$(CC) InvertY.cpp -I $(RINC) -c -o InvertY.o

# PolyaGammaApproxAlt.o: PolyaGammaApproxAlt.h PolyaGammaApproxAlt.cpp
# 	$(CC) PolyaGammaApproxAlt.cpp -I $(RINC) -L $(RLIB) -c -o PolyaGammaApproxAlt.o

# PolyaGamma.o: PolyaGamma.h PolyaGamma.cpp
# 	$(CC) PolyaGamma.cpp -I $(RINC) -L $(RLIB) -c -o PolyaGamma.o

# truncated_norm.o: truncated_norm.h truncated_norm.cpp
# 	$(CC) truncated_norm.cpp -I $(RINC) -L $(RLIB) -c -o truncated_norm.o

# truncated_gamma.o: truncated_gamma.h truncated_gamma.cpp
# 	$(CC) truncated_gamma.cpp -I $(RINC) -L $(RLIB) -c -o truncated_gamma.o

# inverse_gaussian.o: inverse_gaussian.h inverse_gaussian.cpp
# 	$(CC) inverse_gaussian.cpp -I $(RINC) -L $(RLIB) -c -o inverse_gaussian.o

clean:
	rm *.o
