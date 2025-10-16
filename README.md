# **Source code for "Scaling Biclique Percolation for Community Detection in Large Bipartite Networks"**

## **Compile** ##
```
cd bicpc/
cmake .
make
```

## **run MBAG** ##

```
cd bicpc/
bin/BCPC -f ../datasets/youtube.txt -p 2 -q 2 -a baseline -o 0 /* -f filename -p alpha -q beta -a algorithm -o output */
```

## **run PBCPC** ##
```
cd bicpc/
bin/BCPC -f ../datasets/youtube.txt -p 2 -q 2 -a qbicpcLpoor -o 0
```

## **run PBCPC+** ##
```
cd bicpc/
bin/BCPC -f ../datasets/youtube.txt -p 2 -q 2 -a qbicpcL -o 0
```

## **run Biclique** ##
```
cd bicpc/
bin/BCPC -f ../datasets/youtube.txt -p 2 -q 2 -a qpbcl -o 0
```

## **run BicliqueM** ##
```
cd bicpc/
bin/BCPC -f ../datasets/youtube.txt -p 2 -q 2 -a qpbcl_mbc -o 0
```

## **run BicliqueP** ##
```
cd bicpc/
bin/BCPC -f ../datasets/youtube.txt -p 2 -q 2 -a qpbcl_qcpcL -o 0
```