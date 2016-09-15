### Mail from Benjamin -- dated Sep 15 206 ##
## messages 
 
 # Comparision Team,
I finally redid the figures with the corrected data for SC4. I divide the analysis in four parts: number of variants, depth of variants, quality of variants and discordant variants. Below are some results I see from these new figures.

# NUMBERS
In SC1 we find more variants in WGS compared to WES (numbersSC1.pdf). In SC4 the number of variants found is similar between WGS and WES (numbersSC4.pdf). Interestingly when increasing the depth of the WES we do not find more variants exclusive to WES but rather the WES data agrees more with the WGS data (i.e. the ‘Exome only’ bar does not increase from numbersSC1.pdf to numbersSC4.pdf).

The variants found exclusively in WES have many REF variants (numbersRHH2_SC1.pdf and numbersRHH2_SC2.pdf).

# DEPTH
In SC1 discordant variants in WES have low coverage (depthSC1.pdf). For WES in SC4 there is a high variance in the coverage of discordant variants and variants unique to WES have a slightly lower coverage (depthSC4.pdf). 

For WES both HET and HOM variants have higher coverage than REF variants (depthRHH_SC1.pdf and depthRHH_SC4.pdf). HOM variants found exclusively either in WES or WGS have lower coverage (depthRHH2_SC1.pdf and depthRHH2_SC4.pdf).

# QUALITY
For SC1 the quality of discordant variants is low (qualitySC1.pdf). Variants exclusive to either WES or WGS is lower (qualitySC1.pdf and qualitySC4.pdf). 

Quality of HET variants is higher than for other variants (qualityRHH_SC1.pdf and qualityRHH_SC4.pdf). There is more consistency in quality of the WES in SC4 (qualityRHH_SC4.pdf)

There is quite some difference in the trend of the quality for WGS between SC1 and SC4. Not sure why this is - could just be related to the samples themselves.

# DISCORDANT
The majority of discordant variant are the ones called as HET in WGS but REF in WES (discV_SC1.pdf and discV_SC4.pdf). The number of discordant variants are fewer for SC4 than for SC1 (discV_SC1.pdf and discV_SC4.pdf (and numbersSC1.pdf and numbersSC4.pdf)).

Judging from the pattern of the depth and quality of the discordant variants in can be argued that they can be judged as follows:

REF/HET - likely HET as judged by WES due to high coverage of WES data but more so because of low quality of WGS data.

HET/REF ikely HET as judged by WGS due to low quality of WGS data. 

HET/HOM likely HET as judged by WGS due to low quality of WGS data.

HOM/HET is hard to say. Similar quality metrics for both WES and WGS. (But based on the judgments above it seems that HET is more often the correct call.)

REF/HOM untrustworthy due to poor coverage and poor quality in both WES and WGS.

HOM/REF untrustworthy due to poor coverage and poor quality in both WES and WGS.


# HIGHLIGHTS
*30X WGS data is better than 30X WES data (more variants, higher quality).
*Increasing the coverage of the WES data from ~30X to about ~70X improves the WES data but not to the effect that it surpasses the 30X WGS data.
*The depth of HET and HOM variants is much higher than REF variants in WES data.
*The quality of discordant variants is low for WES data.
*The quality of HET variants is higher than for other variants in both WES and WGS.
*The most abundant discordant variants are the ones called as HET in WGS and REF in WES and they can be better trusted as correct in WGS data.
