### Memory and CPU considerations
As per Parse [guidelines](https://support.parsebiosciences.com/hc/en-us/articles/23060102930580-Pipeline-Installation-Current-Version), below is the recommended minimum memory & CPU requirements.  

*For processing a single sublibrary:*

Memory:

- 100M reads or less =  64GB
- 100M-500M reads = 128GB
- 500M-1B  reads = 256GB

CPU/threads:
- 100M reads or less =  8 threads
- 100M-500M reads = 16 threads
- 500M-1B  reads = 24-32 threads

We have 8 sublibraries; each with 1,600M reads after concatenation, so would need aprox 400GB memory & 44 threads


```
cd /exports/eddie/scratch/pdewari/newvolume/

screen -S spipe


```
