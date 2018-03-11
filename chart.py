import sys
from suds import *
from suds.client import Client

def create_chart (genes, chart_path):

    url = 'http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'

    client = Client(url)
    client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')
    #authenticate user email 
    client.service.authenticate('wdong@wdong.org')

    inputIds=','.join(list(genes)) #'201037_at,212561_at,209306_s_at,226137_at,208637_x_at,208862_s_at,200934_at,213457_at,223492_s_at,202651_at'
    categories='GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE'

    idType = 'AFFYMETRIX_3PRIME_IVT_ID'
    listName = 'make_up'
    listType = 0
    client.service.addList(inputIds, idType, listName, listType)
    client.service.setCategories(categories)

    thd = 0.1
    ct = 2
    report = client.service.getChartReport(thd, ct) #overlap, initialSeed, finalSeed, linkage, kappa)

    #parse and print chartReport
    with open(chart_path, 'w') as fOut:
        fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
        for simpleChartRecord in report:
            categoryName = simpleChartRecord.categoryName
            termName = simpleChartRecord.termName
            listHits = simpleChartRecord.listHits
            percent = simpleChartRecord.percent
            ease = simpleChartRecord.ease
            Genes = simpleChartRecord.geneIds
            listTotals = simpleChartRecord.listTotals
            popHits = simpleChartRecord.popHits
            popTotals = simpleChartRecord.popTotals
            foldEnrichment = simpleChartRecord.foldEnrichment
            bonferroni = simpleChartRecord.bonferroni
            benjamini = simpleChartRecord.benjamini
            FDR = simpleChartRecord.afdr
            rowList = [categoryName,termName,str(listHits),str(percent),str(ease),Genes,str(listTotals),str(popHits),str(popTotals),str(foldEnrichment),str(bonferroni),str(benjamini),str(FDR)]
            fOut.write('\t'.join(rowList)+'\n')
            pass
        pass



