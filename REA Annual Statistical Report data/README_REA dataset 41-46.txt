This version includes data in 1938-1946 Rural Electrification Administration (REA) Annual Statistical Reports (ASRs).

I scrape the 1938-1946 REA Annual Statistical Reports to obatin this dataset. It combines variables from several sources：

(1) basic geographic information:
      This includes REA borrower's location variable, like which state, county, city they locate and their coordinate (latitude and longitude). 
      The latitude and longitude of REA site is obtained by two ways (there is a variable called flag_84 to distinguish this)：
      a. From 1984 Active REA Borrower List, scrape the zipcode of REA site, and find latitude and longitude using zipcode;
      b. If the REA site cannot be found in 1984 list, I scrape the town (city) name of REA site in ASRs, and find latitude and longitude using town name.
      Note that the coordinates of REA borrower may be not accurate, since I have not obtained 1940s digital map yet. The coordinates are inferred from 2019 dataset.
      If some zipcode does not exist today, I manually search those coordinates.

(2) Distance to nearest roadshow stop and large city
      We manually collect the location and time (may be accurate to day-level) of roadshow stop from REA news in 1937-1941. Then, we use nncross package of R
      to calculate the distance of each REA borrower to the nearest roadshow stop in each year. Note that since we don't have roadshow stop information from 1942,
      the distance will be the same in all years in 1942-1946. We also calculate the distance of each borrower to the city whose population is larger than 50,000 in 1940.
      City's population data is from 1940 census Databook.

(3) Variable from 1941-1946: 
      This part contains the variables reported year-by-year in ASRs. It includes borrower-year level information including allotment, operation conditon,
      revenue and expenses, and interest payment to government for each REA borrower.  Note that some variables may only be available in specific year,
      and actually the reported variables have undergone a huge change from 1945. Also, from 1945, REA separately report data for distribution and generation
      systems. For more details about each variable, please refer to column "Note" in codebook sheet of dataset.
 
(4) Variable from 1938-1946：
     This part contains the variables reported multiple years in one ASR so that we can data these variables back to 1938.
     This includes miles energized, KWH per residential customer (on December of each year), consumers per mile and revenue per mile (on December of each year).
     Note that these variables are unavailable for generation systems from 1945. 
     To get the yearly level energy consumption data, we merge dataset (3) and (4), and adopt below strategy:
     KWH per residential customer is measured by (KWH billed / Number of Customer) from 1941-1946, 
     while from 1938-1940 is proxied by (KWH per residential customer on December*12) in original REA ASRs report. 
     Also, for 1941-1946, we replace missing values with (KWH per residential customer*12) if KWH billed or number of customer is missing. 

(5) Census and Political Data:
     This part is from ICPSR dataset "Historcal, Demographic, Economics, Social Data, 1790-2002", "United States Agriculture Data, 1840 - 2012"
     and "State Party Strength in the United States: 1872-1996". We use this data as covariates in our propensity score matching analysis. For more datails,
     please refer to their reference manual and codebook.


Reference
[1] US Rural Electrification Administration, Annual Statistical Report 1941-1946, Scanned version of Ohio State University and downloaded from HathiTrust.
[2] Haines, Michael R., and Inter-university Consortium for Political and Social Research. Historical, Demographic, Economic, and Social Data: The United States, 1790-2002. Inter-university Consortium for Political and Social Research [distributor], 2010-05-21. https://doi.org/10.3886/ICPSR02896.v3
[3] Haines, Michael, Fishback, Price, and Rhode, Paul. United States Agriculture Data, 1840 - 2012. Inter-university Consortium for Political and Social Research [distributor], 2018-08-20. https://doi.org/10.3886/ICPSR35206.v4
[4] David, Paul T., and Claggett, William. Party Strength in the United States: 1872-1996. [distributor], 2008-09-10. https://doi.org/10.3886/ICPSR06895.v1