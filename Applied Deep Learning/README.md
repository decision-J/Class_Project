# Personality Detection Method Applied to Korean in Text
** Applied Deep Learning class(prof. Hayoung Kim) Term project



​	In this code, personality detection method is tested whether can be applied to Korean text. As the dataset, title of a Korean article crawled on the Internet is used. Unfortunately, the five personality traits presented above do not match in the news article, so I want to judge whether there is a simple positive trait or not. So, classification is conducted to determine whether the title of article has positive personality, and for this purpose, a sports article that is easy to understand the positive tone. 

​	Although the same architecture as the previous research is used, the absence of data caused by language differences such as emotion lexicon is not reflected. The latter part of the report includes manual tuning to improve model performance.

----

- Dataset(newstitle_기아타이거즈.csv)

Titles of Sports Article : The titles of sports article was collected using crawling. I collected the data at the Naver Sports News tab, the keyword was *KIA Tigers* which is KBO baseball team, and period was from April to September  in 2019 that began last year’s season. It was set to fetch up to 100 news article titles per day, and an average of 550 article titles were collected each month. 



- Python code

**[Term project]_Data_Gain_code** : You can make the dataset by crawling on naver sports news. 

**[Term project]_Source_code** : All functions which are used in the analysis code.

**[Term project]_Analysis_code** : You can conduct the analysis and get results from this code.

 



