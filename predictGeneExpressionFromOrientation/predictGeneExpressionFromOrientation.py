import click
import pandas as pd
import numpy as np
import pybedtools

def getFirstConsecutiveTrues(tfList, reverseList=False):
    tfList=list(tfList)
    if reverseList:
        tfList.reverse()
    outList=[False]*len(tfList)
    for ind,tf in enumerate(tfList):
        if tf:
            outList[ind]=True
        if tf==False:
            break
    if reverseList:
        outList.reverse()
    return outList
def getNearbyGeneAnnots(IDstring,geneBedFile,chromsizesFile,slopSize=5000):
    chrom=IDstring.split(':')[0]
    start=int(IDstring.split(':')[1].split('-')[1])
    end=int(IDstring.split(':')[1].split('-')[2])

    pgapBed=pybedtools.bedtool.BedTool(geneBedFile)
    testBed=pybedtools.bedtool.BedTool(f"{chrom} {start} {end}",from_string=True)
    df=pgapBed.intersect(testBed.slop(b=slopSize,g=chromsizesFile),wa=True).to_dataframe()
    if len(df)>0:
        df['ID']=IDstring
        df['upstream']=(df.end<start)&(df.strand=="-")
        df['downstream']=(df.start>end)&(df.strand=="+")
        df['intragenic']=(df.start<=start)&(df.end>=end)
        df['partialstart']=(~df.intragenic)&(((df.end>=start)&(df.start<=start)&(df.strand=="-"))|((df.end>=end)&(df.start<=end)&(df.strand=="+")))
        df['partialstop']=(~df.intragenic)&(((df.end>=start)&(df.start<=start)&(df.strand=="+"))|((df.end>=end)&(df.start<=end)&(df.strand=="-")))
        if (df.end<start).sum():
            df.loc[df.end<start,'upstream']=getFirstConsecutiveTrues(df.upstream[df.end<start],reverseList=True)
        if (df.start<end).sum():
            df.loc[df.start>end,'downstream']=getFirstConsecutiveTrues(df.downstream[df.start>end],reverseList=False)

        df['type']='other'
        df.loc[df.upstream|df.downstream,'type']='instream'
        df.loc[df.intragenic|df.partialstart|df.partialstop,'type']='genic'
    if len(df)>0:
        return df
    else:
        return pd.DataFrame([])

@click.command()
@click.option('-i','--invertons', required=True, help='file listing inverton IDs',type=click.Path())
@click.option('-g','--genes', required=True, help='bed file listing gene locations and strand information',type=click.Path())
@click.option('-p','--promoters', required=True, help='bed file listing promoter locations and strand information',type=click.Path())
@click.option('-c','--chromsizes', required=True, help='chromsizes file',type=click.Path())
@click.option('-s','--slop', default=5000, help='how far in bp to look for genes near invertons (bedtools slop parameter), default 5000bp',type=int)
@click.option('-o','--output', required=True, help='output file name',type=click.Path())
def predict(invertons, genes, promoters, chromsizes, slop, output):
    invertonIDDf=pd.read_csv(invertons,header=None,names=['ID'])
    invertonIDDf['chrom']=[IDstring.split(':')[0] for IDstring in invertonIDDf.ID]
    invertonIDDf['start']=[int(IDstring.split(':')[-1].split('-')[0]) for IDstring in invertonIDDf.ID]
    invertonIDDf['end']=[int(IDstring.split(':')[-1].split('-')[-1]) for IDstring in invertonIDDf.ID]
    invertonBed=pybedtools.BedTool.from_dataframe(invertonIDDf[['chrom','start','end','ID']].sort_values(by=['chrom','start','end']))
    nearbyGeneDf=pd.concat([getNearbyGeneAnnots(ID,geneBedFile=genes,chromsizesFile=chromsizes,slopSize=slop) for ID in invertonIDDf.ID])
    promoterInvertonDf=pybedtools.BedTool(promoters).intersect(invertonBed,wo=True,f=1).to_dataframe()
    promoterInvertonDf=promoterInvertonDf[['chrom','start','end','name','score','strand','blockCount']].rename(columns={'blockCount':'ID'})
    geneExpressionPredictionDf=promoterInvertonDf.merge(nearbyGeneDf.loc[nearbyGeneDf['type']=='instream',['start','end','name','score','strand','ID']],on='ID',suffixes=['_promoter','_gene'])
    geneExpressionPredictionDf['expressedOrientation']=['F' if strand_promoter==strand_gene else 'R' for strand_promoter,strand_gene in zip(geneExpressionPredictionDf.strand_promoter,geneExpressionPredictionDf.strand_gene)]
    geneExpressionPredictionDf=geneExpressionPredictionDf[['chrom','start_promoter','end_promoter','name_promoter','score_promoter','strand_promoter','start_gene','end_gene','name_gene','score_gene','strand_gene','ID','expressedOrientation']]
    geneExpressionPredictionDf.columns=['chrom','start_promoter','end_promoter','name_promoter','score_promoter','strand_promoter','start_gene','end_gene','name_gene','score_gene','strand_gene','inverton ID','expressedOrientation']
    geneExpressionPredictionDf.to_csv(output,sep="\t",index=False)
if __name__ == '__main__':
    predict()