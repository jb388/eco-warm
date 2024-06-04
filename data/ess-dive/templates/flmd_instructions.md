# Instructions for ESS-DIVE File-level Metadata Reporting Format

1. Create file-level metadata for each data set using the provided template  
&nbsp;&nbsp;&nbsp;a. Start new row for each file and include all files  
&nbsp;&nbsp;&nbsp;b. Use "\*" wildcard when the FLMD applies to multiple files  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;i. For example - the same FLMD applies to all soil core files in this data package - "soil\_cores\_\*\.csv"  
2. Follow the instructions in the [FLMD quick guide](flmd_quick_guide.md)  
3. Use the [FLMD template](https://github.com/ess-dive-community/essdive-file-level-metadata/blob/master/flmd_template.xlsx)  
4. Save the FLMD template as a CSV following the CSV Reporting Format guidance. We recommend naming the flmd file "flmd.csv" or you can include a prefix in the form of "\*\_flmd.csv"  
5. For CSV data files, create a [CSV Data Dictionary](CSV_dd/) to describe the fields and other attributes of your CSV data file. Also list any data dictionaries you create as part of your file-level metadata (described in step 1a) 
6. For zip files, we recommend including the FLMD outside the zip file but may also include it in the zip file  

**Notes**

* Fields common to ESS-DIVE Package Level Metadata are consistent in format structure  

**Contents of the Elements**

* File\_Name  
* File\_Description  
* Standard  
* UTC\_Offset  
* File\_Version  
* Contact  
* Date\_Start  
* Date\_End  
* Northwest\_Latitude\_Coordinate  
* Northwest\_Longitude\_Coordinate  
* Southeast\_Latitude\_Coordinate  
* Southeast\_Longitude\_Coordinate  
* Latitude  
* Longitude  
* Missing\_Value\_Codes    
* Data_Orientation  
* Notes

