rm results.txt 
touch results.txt 
rm cuda.txt
touch cuda.txt
rm kernel.txt
touch kernel.txt

nvcc -O3 CUDA_Any_Char_Range_Count_std.cu -o CUDA_Any_Char_Range_Count_std.out

from=a

exe="./CUDA_Any_Char_Range_Count_std.out"
n=( 65536 131072 262144 524288 1048576 2097152 4194304 8388608 16777216 33554432 67108864 134217728 268435456 )

	for runs in {1..12}
		do
		for ii in "${n[@]}"
			 #echo "n=" $n "\n"
			 do
			 k=1
			 y=a
			 for k in {1..26}
				  do
				  #perf stat -e instructions -e cycles -e branch-instructions -e branch-misses -e cache-references -e cache-misses $exe $n $from $y >> results.txt
				$exe $ii $from $y >> results.txt
				  y=$(echo "$y" | tr "0-9a-z" "1-9a-z_")
				  k=$((1 + $k))
				  done
			 
			 done
	done

grep "Time" results.txt | cut -d ':' -f2 | tr -s ' ' > CPU.txt
grep "Time CUDA" results.txt | cut -d ':' -f2 | tr -s ' ' > cuda.txt
grep "kernel ONLY" results.txt | cut -d ':' -f2 | tr -s ' ' > kernel.txt
#grep "Time pthreads" results.txt | cut -d ':' -f2 | tr -s ' ' > phreads.txt
#grep "Time pthreads dynamic" results.txt | cut -d ':' -f2 | tr -s ' ' > phreadsdyn.txt
#grep "Time pthreads AVX2" results.txt | cut -d ':' -f2 | tr -s ' ' > phreadsAVX2.txt
#grep "Time pthreads AVX2 dynamic" results.txt | cut -d ':' -f2 | tr -s ' ' > phreadsAVX2dyn.txt
#grep "Time omp CPU" results.txt | cut -d ':' -f2 | tr -s ' ' > ompCPU.txt
#grep "Time omp CPU dynamic" results.txt | cut -d ':' -f2 | tr -s ' ' > ompCPUdyn.txt
#grep "Time omp AVX2" results.txt | cut -d ':' -f2 | tr -s ' ' > ompAVX2.txt
#grep "Time omp AVX2 dynamic" results.txt | cut -d ':' -f2 | tr -s ' ' > ompAVX2dyn.txt
#grep "Time AVX2" results.txt | cut -d ':' -f2 | tr -s ' ' > avx.txt
#grep ' instructions' text.txt | cut -d '#' -f1 | tr -d 'instructions' | tr -d ',' | tr -s ' '  > instructions.txt
#grep 'branch-instructions' text.txt | cut -d '#' -f1 | tr -d 'branch-instructions' | tr -d ',' | tr -s ' ' | tr -d '-' > br-instructions.txt
#grep 'cycles' text.txt | cut -d '#' -f1 | tr -d 'cycles' | tr -d ',' | tr -s ' ' > cycles.txt
#grep 'cache-references' text.txt | cut -d '' -f1 | tr -d 'cache-references' | tr -d ',' | tr -s ' ' | tr -d '-' > cache-ref.txt
#grep 'cache-misses' text.txt | cut -d '#' -f1 | tr -d 'cache-misses' | tr -d ',' | tr -s ' ' | tr -d '-' > cache-misses.txt
#grep 'branch-misses' text.txt | cut -d '#' -f1 | tr -d 'branch-misses' | tr -d ',' | tr -s ' ' | tr -d '-' > br-misses.txt


