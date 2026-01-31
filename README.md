# Metagenomic Assemblers Comparison
Klára Kostelníková

# Úvod

Cílem této práce je porovnat výkon dvou moderních metagenomických
assemblerů: **metaMDBG** a **myloasm**. Porovnání je prováděno na
reálných datech z Nanopore sekvenování mikrobiálního společenstva z
horkého pramene (hot spring).

Analýza se zaměřuje na vyhodnocení kvality sestavených genomů bez
nutnosti spouštět samotné skládání nebo mapování. Využíváme textová
metadata získaná z:

1.  **Hlaviček contigů** (délka, pokrytí, cirkularita).
2.  **CheckM2** (odhad kompletnosti a kontaminace).
3.  **GTDBtk** (taxonomická klasifikace).

Celé zpracování dat, vizualizace výsledků a finální vyhodnocení je
provedeno v prostředí **R** s využitím balíčků **tidyverse**. Hlavním
kritériem úspěšnosti je schopnost assemblerů rekonstruovat velké (\> 500
kbp) a vysoce kvalitní cirkulární contigy.

# 1. Příprava dat

V této části načítám potřebné knihovny a data.

``` r
library(tidyverse)
theme_set(theme_minimal())
```

## Načtení a čištění dat

### Myloasm

``` r
# 1. Načteme surový text
raw_myloasm <- read_lines("results/myloasm_assembly_headers.txt")

# 2. Vytáhneme data
headers_myloasm <- tibble(raw = raw_myloasm) %>%
  extract(
    col = raw,
    into = c("contig_id", "length", "circ_status", "coverage"),
    # Regex hledá vzory jako "_len-" nebo "_depth-"
    regex = ">([^_]+)_len-([0-9]+)_circular-([a-z]+)_depth-([0-9.]+)",
    convert = TRUE
  ) %>%
  mutate(
    assembler = "myloasm",
    #"possible" = circular
    is_circular = if_else(circ_status %in% c("yes", "possible"), TRUE, FALSE)
  ) %>%
  select(contig_id, assembler, length, coverage, is_circular)

# KONTROLA
head(headers_myloasm)
```

    # A tibble: 6 × 5
      contig_id   assembler length coverage is_circular
      <chr>       <chr>      <int>    <int> <lgl>      
    1 u3840050ctg myloasm    10849        2 FALSE      
    2 u418451ctg  myloasm    23116        2 FALSE      
    3 u911112ctg  myloasm    78285       27 FALSE      
    4 u2187011ctg myloasm    16454        2 FALSE      
    5 u3322747ctg myloasm    14907        1 FALSE      
    6 u2115919ctg myloasm     9865        1 FALSE      

### MetaMDBG

``` r
# 1. Načteme surový text
raw_metamdbg <- read_lines("results/metamdbg_assembly_headers.txt")

# 2. Vytáhneme data 
headers_metamdbg <- tibble(raw = raw_metamdbg) %>%
  mutate(
    # Najde ID (první slovo po >)
    contig_id = str_extract(raw, "(?<=^>)\\S+"),
    # Najde číslo za "len="
    length = str_extract(raw, "(?<=length=)[0-9]+") %>% as.numeric(),
    # Najde číslo za "depth="
    coverage = str_extract(raw, "(?<=coverage=)[0-9.]+") %>% as.numeric(),
    # Pokud tam je "circular=yes", dá TRUE
    is_circular = str_detect(raw, "circular=yes"),
    assembler = "metamdbg"
  ) %>%
  select(contig_id, assembler, length, coverage, is_circular)

# KONTROLA
head(headers_metamdbg)
```

    # A tibble: 6 × 5
      contig_id assembler length coverage is_circular
      <chr>     <chr>      <dbl>    <dbl> <lgl>      
    1 ctg0      metamdbg    7710     2.25 FALSE      
    2 ctg1      metamdbg    8229     1.96 FALSE      
    3 ctg2      metamdbg   10850     3.32 FALSE      
    4 ctg3      metamdbg    9023     4.85 FALSE      
    5 ctg4      metamdbg   56052     7.9  FALSE      
    6 ctg5      metamdbg   12495     2.07 FALSE      

### Načtení kvality

``` r
# 1. Spojíme hlavičky pod sebe
headers_all <- bind_rows(headers_myloasm, headers_metamdbg)

# 2. Načteme CheckM2 reporty (kvalitu)
qual_myloasm <- read_tsv("results/checkm2/myloasm/quality_report_myloasm.tsv", show_col_types = FALSE) %>%
  select(Name, Completeness, Contamination) %>%
  rename(contig_id = Name) 

qual_metamdbg <- read_tsv("results/checkm2/metamdbg/quality_report.tsv", show_col_types = FALSE) %>%
  select(Name, Completeness, Contamination) %>%
  rename(contig_id = Name)

quality_all <- bind_rows(qual_myloasm, qual_metamdbg)

#KONTROLA:
head(quality_all)
```

    # A tibble: 6 × 3
      contig_id   Completeness Contamination
      <chr>              <dbl>         <dbl>
    1 u1005671ctg         50.8          0.6 
    2 u1022233ctg         95.4          0.4 
    3 u1027066ctg         89.9          2.92
    4 u1041087ctg         11.3          0.02
    5 u1047122ctg         91.0          0.21
    6 u1079443ctg         95.6          0.18

Načtení taxonomie z GTDBtk

``` r
# 1. Načtení GTDBtk taxonomie

tax_bac_ <- read_tsv("results/gtdbtk/myloasm/classify/gtdbtk.bac120.summary.tsv", show_col_types = FALSE) %>% 
  select(user_genome, classification)

tax_bac_metamdbg <- read_tsv("results/gtdbtk/metamdbg/classify/gtdbtk.bac120.summary.tsv", show_col_types = FALSE) %>% 
  select(user_genome, classification)

tax_bac <- bind_rows(tax_bac_metamdbg, tax_bac_)
  
tax_ar_ <- read_tsv("results/gtdbtk/myloasm/classify/gtdbtk.ar53.summary.tsv", show_col_types = FALSE) %>% 
  select(user_genome, classification)

tax_ar_metamdbg <- read_tsv("results/gtdbtk/metamdbg/classify/gtdbtk.ar53.summary.tsv", show_col_types = FALSE) %>% 
  select(user_genome, classification)

tax_ar <- bind_rows(tax_ar_metamdbg, tax_ar_)

# 2. Spojíme bakterie a archea do jedné tabulky a rozdělíme klasifikaci
taxonomy_all <- bind_rows(tax_bac, tax_ar) %>%
  rename(contig_id = user_genome) %>%
  # Rozdělíme sloupec 'classification' na 7 nových sloupců podle středníku
  separate(
    classification, 
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"), 
    sep = ";", 
    fill = "right"
  ) %>%
  # Odstraníme prefixy (např. "p__", "c__") ze všech nově vytvořených sloupců
  mutate(across(c(domain:species), ~ str_remove(., "^[a-z]__"))) %>%
  # Odstranění řádků, kde chybí phylum 
  filter(!is.na(phylum) & phylum != "")
```

spojení dat

``` r
# 1.Čištění hlaviček
headers_clean <- headers_all %>%
  mutate(contig_id = str_replace(contig_id, "\\.fasta.*", "")) %>%  
  mutate(contig_id = if_else(
    str_starts(contig_id, "u") & str_detect(contig_id, "_"),
    str_extract(contig_id, "^[^_]+"), 
    contig_id)) %>%
  mutate(contig_id = str_trim(contig_id))

# 2. Očištění kvality (CheckM2)
quality_clean <- quality_all %>%
  mutate(contig_id = str_replace(contig_id, "\\.fasta.*", "")) %>%  
  mutate(contig_id = if_else(
    str_starts(contig_id, "u") & str_detect(contig_id, "_"),
    str_extract(contig_id, "^[^_]+"), 
    contig_id )) %>%
  mutate(contig_id = str_trim(contig_id))

taxonomy_clean <- taxonomy_all %>%
  mutate(contig_id = str_replace(contig_id, "\\.fasta.*", "")) %>%  
  mutate(contig_id = if_else(
    str_starts(contig_id, "u") & str_detect(contig_id, "_"),
    str_extract(contig_id, "^[^_]+"), 
    contig_id )) %>%
  mutate(contig_id = str_trim(contig_id))



  


# --- FINÁLNÍ SPOJENÍ ---

final_data <- headers_clean %>%
  # Připojíme kvalitu z CheckM2
  inner_join(quality_clean, by = "contig_id") %>%
  # Připojíme taxonomii z GTDBtk 
  inner_join(taxonomy_clean, by = "contig_id")

# Kontrola výsledku
glimpse(final_data)
```

    Rows: 130
    Columns: 14
    $ contig_id     <chr> "u1859167ctg", "u3662534ctg", "u3457852ctg", "u1614917ct…
    $ assembler     <chr> "myloasm", "myloasm", "myloasm", "myloasm", "myloasm", "…
    $ length        <dbl> 2142282, 1109361, 1142602, 1752029, 1032745, 519055, 652…
    $ coverage      <dbl> 10, 5, 164, 30, 6, 6, 7, 7, 6, 7, 11, 6, 7, 7, 10, 7, 6,…
    $ is_circular   <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, …
    $ Completeness  <dbl> 99.67, 43.81, 96.70, 81.74, 99.97, 15.23, 56.08, 50.78, …
    $ Contamination <dbl> 0.66, 0.07, 0.00, 4.75, 0.22, 0.01, 0.01, 0.60, 0.00, 1.…
    $ domain        <chr> "Bacteria", "Bacteria", "Archaea", "Bacteria", "Bacteria…
    $ phylum        <chr> "UBA6262", "Nitrospirota", "Nanobdellota", "Firestonebac…
    $ class         <chr> "UBA6262", "Thermodesulfovibrionia", "Nanobdellia", "", …
    $ order         <chr> "UBA6262", "Thermodesulfovibrionales", "Woesearchaeales"…
    $ family        <chr> "UBA6262", "SM23-35", "JACPBV01", "", "GCA-002772895", "…
    $ genus         <chr> "UBA6262", "JBFMAM01", "", "", "", "JACPQH01", "JBFLZU01…
    $ species       <chr> "", "", "", "", "", "", "JBFLZU01 sp040756155", "", "", …

# 2. Výsledky a Vizualizace

## Distribuce délky contigů

``` r
final_data <- final_data %>%
  mutate(quality_cat = case_when(
    Completeness > 90 & Contamination < 5 ~ "High",
    Completeness > 50 & Contamination < 10 ~ "Medium",
    TRUE ~ "Low"
  ))

# Rychlé ověření, kolik jich máme v každé kategorii
table(final_data$assembler, final_data$quality_cat)
```

              
               High Low Medium
      metamdbg   34   6      2
      myloasm    63  14     11

**Korelace délky a pokrytí:** Z bodového grafu (scatter plot) je patrné,
že většina cirkulárních contigů se nachází v oblasti středního až
vysokého pokrytí. U assembleru **myloasm** vidíme, že dokázal
zrekonstruovat dlouhé contigy i při nižším pokrytí, zatímco **metaMDBG**
vyžadoval pro stabilní sestavení dlouhých sekvencí spíše vyšší hloubku
čtení. Lineární trend zde není striktní, což naznačuje, že délka
sestaveného genomu (MAG) v tomto vzorku nezávisí pouze na množství dat,
ale i na komplexitě daného organismu

``` r
final_data %>%
  ggplot(aes(x = length, fill = is_circular)) +
  geom_histogram(bins = 50, color = "black", alpha = 0.7) +
  scale_x_log10(labels = scales::comma) + # Logaritmická osa X
  facet_wrap(~assembler, ncol = 1) +
  labs(title = "Distribution of contig lenght",
       x = "Lenght (bp, log scale)", y = "Number", fill = "Circular")
```

![](worflow_files/figure-commonmark/plot-length-1.png)

``` r
# Boxplot pro srovnání délek
final_data %>%
  ggplot(aes(x = assembler, y = length, fill = is_circular)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "Contig Length Distribution by Assembler",
       subtitle = "Comparison of circular vs. non-circular contigs (log scale)",
       x = "Assembler",
       y = "Contig Length (bp)",
       fill = "Circular") +
  theme_minimal()
```

![](worflow_files/figure-commonmark/plot-length-2.png)

## Kvalita MAGs

Zajímá nás kvalita velkých cirkulárních contigů.

``` r
# Vytvoření výběru 
large_circular <- final_data %>%
  filter(length > 500000, is_circular == TRUE)


# Number of MAGs per phylum colored by quality category
large_circular %>%
  filter(!is.na(phylum)) %>% 
  count(assembler, phylum, quality_cat) %>% 
  ggplot(aes(x = phylum, y = n, fill = quality_cat)) + 
  geom_col() +
  coord_flip() + 
  facet_wrap(~assembler) +
  scale_fill_manual(values = c("High" = "#2ca25f", "Medium" = "#feb24c", "Low" = "#de2d26")) +
  labs(title = "Number of Large Circular MAGs per Phylum",
       subtitle = "Contigs > 500 kb categorized by quality",
       x = "Phylum", 
       y = "Count", 
       fill = "Quality") +
  theme_minimal()
```

![](worflow_files/figure-commonmark/plot-quality-1.png)

``` r
# Bodový graf kvality s referenčními čarami a popisky
large_circular %>%
  ggplot(aes(x = Contamination, y = Completeness, color = quality_cat)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~assembler) +
  # Referenční čáry pro High Quality (90% completeness, 5% contamination)
  geom_hline(yintercept = 90, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey40") +
  # Nastavení barev pro kategorie
  scale_color_manual(values = c("High" = "#2ca25f", "Medium" = "#feb24c", "Low" = "#de2d26")) +
 labs(title = "Completeness and Contamination of Large MAGs",
       subtitle = "Circular contigs > 500 kb; Dashed lines indicate High Quality thresholds",
       x = "Contamination (%)", 
       y = "Completeness (%)",
       color = "Quality Category") +
  theme_minimal()
```

![](worflow_files/figure-commonmark/plot-quality-2.png)

# Závěr

V této práci jsme porovnali výkon assemblerů **metaMDBG** a **myloasm**
při zpracování nanopore dat z horkého pramene.

- **Strukturální integrita:** Assembler **myloasm** vykazuje vyšší
  efektivitu v rekonstrukci dlouhých sekvencí, což dokládá distribuce
  délek contigů. Významným rozdílem je schopnost identifikovat “possible
  circular” contigy, které po filtraci nad 500 kbp poskytují ucelenější
  pohled na genomickou informaci.

- **Kvalita MAGs:** Dle kritérií CheckM2 oba nástroje generují primárně
  vysoce kvalitní genomy (\>90 % kompletnost, \<5 % kontaminace).
  Nicméně **myloasm** doručil vyšší absolutní počet těchto vysoce
  kvalitních (High Quality) cirkulárních MAGs.

- **Taxonomická reprezentace:** Analýza pomocí GTDB-Tk odhalila
  dominantní zastoupení kmene **Patescibacteriota**. **Myloasm** úspěšně
  zrekonstruoval širší spektrum zástupců, což naznačuje jeho lepší
  citlivost k méně zastoupeným nebo specifickým taxonomickým skupinám v
  daném prostředí.

  **Celkové zhodnocení:** Na základě počtu a kvality zrekonstruovaných
  velkých cirkulárních contigů lze konstatovat, že pro analýzu tohoto
  mikrobiálního společenstva je vhodnějším nástrojem **myloasm**.
  Poskytuje robustnější výsledky pro následné metabolické a funkční
  analýzy
