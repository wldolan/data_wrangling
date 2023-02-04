library("easyPubMed")
library(tidyverse)
library(readxl)
library(googlesheets4)
library(cronR)

### googlsheets authentication
gs4_auth(
  email = "whitney@perlara.com",
  path = NULL,
  scopes = "https://www.googleapis.com/auth/spreadsheets",
  cache = gargle::gargle_oauth_cache(),
  use_oob = gargle::gargle_oob_default(),
  token = NULL
)


#### 8p query
gs_last <- read_sheet("https://docs.google.com/spreadsheets/d/17F9mAOJiSC_70oElnu1r2-woBRvjJHEpjSdasixO-04/",
                            sheet = "8p", skip = 1)
#last_run <- "2022/01/01"
last_run <- (max(gs_last$Run_id)) %>% str_replace_all("-", "/") 

query <- paste0("(\"", last_run,"\"[Date - Publication] : \"3000\"[Date - Publication]) AND (8p AND (cnv OR rearrangement OR invdupdel OR inv-dup OR inversion OR deletion OR duplication))")
my_query <- get_pubmed_ids(query)

if (my_query$Count == "0") {
  print("no results") 
} else {
  
  # Fetch data
  my_abstracts_xml <- fetch_pubmed_data(my_query, retmax = 5000)
  
  # Store Pubmed Records as elements of a list
  all_xml <- articles_to_list(my_abstracts_xml)
  
  # Perform operation (use lapply here, no further parameters)
  df <- do.call(rbind, lapply(all_xml, article_to_df,
                              max_chars = -1,
                              getAuthors = TRUE,
                              getKeywords = TRUE
  ))
  
  # reshape
  authors_df <- df %>%
    unite("author", c(lastname, firstname), remove=TRUE) %>%
    group_by(pmid) %>%
    summarize(authors = paste0(author, collapse=", ")) 
  
  df_out <- df %>%
    unite("pubdate", c(year, month, day), sep = "-", remove=TRUE) %>%
    select(-c(address, email, lastname, firstname)) %>%
    distinct() %>%
    arrange(pubdate) %>%
    left_join(authors_df, by = "pmid") %>%
    mutate(query = query) %>%
    mutate(run_id = Sys.Date())
  
  # export to googlesheet
  gs <- gs4_get("https://docs.google.com/spreadsheets/d/17F9mAOJiSC_70oElnu1r2-woBRvjJHEpjSdasixO-04/")
  
  sheet_append(gs, df_out, sheet = "8p")
}

#### repeat everything for gene queries ####
gs_last <- read_sheet("https://docs.google.com/spreadsheets/d/17F9mAOJiSC_70oElnu1r2-woBRvjJHEpjSdasixO-04/",
                      sheet = "genes", skip=1)

last_run <- max(gs_last$Run_id) %>% str_replace_all("-", "/") 
#last_run <- "2023/01/01"

input_file <- read_xlsx("~/Documents/Perlara/8P/8p Gene Information w_ Phenotype Scores.xlsx",
                        sheet = "gene_subset", skip =1)

### search and compile complete results
genelist <- as.list(input_file$gene)
genelist[39] <- "HR gene" 

datalist <- list()
n <- 0
system.time(
  for (g in 1:length(genelist)){
    gene <- genelist[g]
    query_part1 <- paste0("(\"", last_run,"\"[Date - Publication] : \"3000\"[Date - Publication]) AND (gene[TIAB] AND (cognitive[TIAB] OR neuro[TIAB] OR pathogenesis[TIAB] OR disease[TIAB] OR cardio[TIAB] OR development[TIAB] OR therapeutic[TIAB])") 
    query_part2 <- paste0(" AND ", gene, ")")
    gene_query <- paste(query_part1, query_part2, sep = " ") 
    
    my_query <- get_pubmed_ids(gene_query)
    
    if (my_query$Count == "0") {
      next
    }
    
    # Fetch data
    my_abstracts_xml <- fetch_pubmed_data(my_query, retmax = 5000)
    
    # Store Pubmed Records as elements of a list
    all_xml <- articles_to_list(my_abstracts_xml)
    
    # Perform operation (use lapply here, no further parameters)
    df <- do.call(rbind, lapply(all_xml, article_to_df,
                                max_chars = -1,
                                getAuthors = TRUE,
                                getKeywords = TRUE
    ))
    
    df$gene <- as.character(gene)
    df$query <- as.character(query)

    n <- n+1
    datalist[[n]] <- df
    
    print(paste('gene', g, gene))
      
  }
)

genes_df <- do.call(rbind, datalist)

# reshape
authors_df <- genes_df %>%
  unite("author", c(lastname, firstname), remove=TRUE) %>%
  group_by(pmid) %>%
  summarize(authors = paste0(author, collapse=", ")) 

df_out <- genes_df %>%
  unite("pubdate", c(year, month, day), sep = "-", remove=TRUE) %>%
  select(-c(address, email, lastname, firstname)) %>%
  distinct() %>%
  arrange(pubdate) %>%
  left_join(authors_df, by = "pmid") %>%
  mutate(run_id = Sys.Date()) %>%
  select(pmid, doi, title, abstract, pubdate, jabbrv, journal, keywords, authors, gene, query, run_id) %>%
  filter(!str_detect(tolower(journal), "cancer|onco|veterinary|plant|fish|bovine|microb|botany|STAR protocols|animal|environment"))

# export to googlesheet
gs <- gs4_get("https://docs.google.com/spreadsheets/d/17F9mAOJiSC_70oElnu1r2-woBRvjJHEpjSdasixO-04/")

sheet_append(gs, df_out, sheet = "genes")
