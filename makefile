xsir_evo_plast_gandon : host_specialization_plasticity.cpp auxiliary.h
	g++ -Wall -O3 -o xhost_specialization_plasticity host_specialization_plasticity.cpp  -lm -lrt -lgsl -lgslcblas

.PHONY: clean

clean:
	rm -rf xhost_specialization_plasticity 
