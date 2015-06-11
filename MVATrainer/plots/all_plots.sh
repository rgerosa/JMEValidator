mkdir plots/$1
for f in plots/plots/input*.json
do
 echo "Processing $f"
 harry.py -j $f -o plots/$1/input
done

for f in plots/plots/corrected*.json
do
 echo "Processing $f"
 harry.py -j $f -o plots/$1/corrected
done

for f in plots/plots/mva*.json
do
 echo "Processing $f"
 harry.py -j $f -o plots/$1/mva
done

for f in plots/plots/performance*.json
do
 echo "Processing $f"
 harry.py -j $f -o plots/$1/performance
done

for f in plots/plots/cov*.json
do
 echo "Processing $f"
 harry.py -j $f -o plots/$1/cov
done
