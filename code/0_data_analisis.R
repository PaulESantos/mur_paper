# Libraries ---------------------------------------------------------------
library(tidyverse)
library(here)
library(janitor)
library(readxl)
library(ggview)
library(vegan)
source("code/funciones.R")
# Importar data ------------------------------------------------------------

data <- readxl::read_excel("data/base de datos final.xlsx", sheet = 1) |>
  janitor::clean_names() |>
  dplyr::select(fecha, especie, tipo_de_bosque) |>
  tidyr::drop_na() |>
  dplyr::mutate(especie = stringr::str_trim(especie) |>
    stringr::str_to_sentence())
data
abrev <- readxl::read_excel("data/base de datos final.xlsx", sheet = 2) |>
  janitor::clean_names() |>
  dplyr::mutate(especie = stringr::str_trim(especie) |>
    stringr:::str_to_sentence()) |>
  dplyr::rename(species = especie)
abrev

# Por tipo de bosque --------------------------------------------

comm_mur <- data |>
  dplyr::select(tipo_de_bosque, especie) |>
  dplyr::mutate(especie = stringr::str_to_sentence(especie)) |>
  dplyr::group_by(tipo_de_bosque, especie) |>
  dplyr::summarise(
    abundancia = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(tipo_de_bosque = dplyr::case_when(
    tipo_de_bosque == "Bosque de colina baja" ~ "BCB",
    tipo_de_bosque == "Bosque de montaña basimontano" ~ "BMBM",
    tipo_de_bosque == "Bosque de montaña montano" ~ "BMM",
    tipo_de_bosque == "Bosque de terraza baja" ~ "BTB"
  ))

comm_mur |>
  pesa::check_na()

comm_mur <- comm_mur |>
  tbl_to_comm(tipo_de_bosque, especie, abundancia)

comm_mur

# Tipo de bosque por dia

comm_dia <- data |>
  dplyr::select(fecha, tipo_de_bosque, especie) |>
  dplyr::mutate(daya = lubridate::day(fecha)) |>
  dplyr::group_by(tipo_de_bosque, daya, especie) |>
  dplyr::summarise(
    abun = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    tipo_de_bosque = dplyr::case_when(
      tipo_de_bosque == "Bosque de colina baja" ~ "BCB",
      tipo_de_bosque == "Bosque de montaña basimontano" ~ "BMBM",
      tipo_de_bosque == "Bosque de montaña montano" ~ "BMM",
      tipo_de_bosque == "Bosque de terraza baja" ~ "BTB"
    ),
    tipo_de_bosque = paste0(tipo_de_bosque, "_", daya)
  )

comm_dia |>
  tbl_to_comm(tipo_de_bosque, especie, abun)

## Bosque de colina baja

bcb_comm <- comm_dia |>
  dplyr::filter(stringr::str_detect(tipo_de_bosque, "BCB")) |>
  tbl_to_comm(tipo_de_bosque, especie, abun)
bcb_comm

## Bosque de montaña basimontano

bmbm_comm <- comm_dia |>
  dplyr::filter(stringr::str_detect(tipo_de_bosque, "BMBM")) |>
  tbl_to_comm(tipo_de_bosque, especie, abun)

## Bosque de montaña montano

bmm_comm <- comm_dia |>
  dplyr::filter(stringr::str_detect(tipo_de_bosque, "BMM")) |>
  tbl_to_comm(tipo_de_bosque, especie, abun)

## Bosque de terraza baja

btb_comm <- comm_dia |>
  dplyr::filter(stringr::str_detect(tipo_de_bosque, "BTB")) |>
  tbl_to_comm(tipo_de_bosque, especie, abun)



# Curva de acumulación ------------------------------------------

accum_sitio <- dplyr::bind_rows(get_accum_data(bcb_comm, method = "exact"),
  get_accum_data(bmbm_comm, method = "exact"),
  get_accum_data(bmm_comm, method = "exact"),
  get_accum_data(btb_comm, method = "exact"),
  .id = c("site")
) |>
  dplyr::mutate(site = dplyr::case_when(
    site == "1" ~ "BCB",
    site == "2" ~ "BMBM",
    site == "3" ~ "BMM",
    site == "4" ~ "BTB"
  ))
accum_sitio

accum_sitio |>
  ggplot(aes(sites, richness, color = site, fill = site)) +
  geom_ribbon(
    aes(
      ymin = max,
      ymax = min,
      fill = site
    ),
    color = "transparent",
    alpha = .2
  ) +
  geom_line(aes(group = site), linetype = "dotted", linewidth = 1.2) +
  geom_linerange(aes(
    ymin = min,
    ymax = max
  )) +
  geom_point(size = 3, shape = 19) +
  scale_y_continuous(
    limits = c(
      min(accum_sitio$min),
      (max(accum_sitio$richness) + max(accum_sitio$sd))
    ),
    breaks = seq(0, 40, 5)
  ) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  theme_bw() +
  labs(
    x = "Lugares",
    y = str_to_sentence(unique(accum_sitio$method)),
    title = paste0(
      "Curva de Acumulación de Especies - ",
      str_to_sentence(unique(accum_sitio$method))
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  labs(
    fill = "Bosque", color = "Bosque",
    x = "Número de días",
    y = "Especies"
  )

ggview(
  dpi = 600,
  units = "cm",
  height = 9,
  width = 15
)
ggsave("img/species_acumulation.png",
  dpi = 600,
  units = "cm",
  height = 9,
  width = 15
)

# Diversidad alfa -----------------------------------------------

alpha_mur <- adiv::speciesdiv(comm_mur) |>
  dplyr::as_tibble(rownames = "sites") |>
  janitor::clean_names() |>
  dplyr::mutate(pielou = shannon / log(richness))

alpha_mur

# Curva rango abundancia ----------------------------------------


rank_abund_data <- rankabund_df(comm_mur, group = "sites") |>
  ungroup()
rank_abund_data

rank_abund_data1 <- rank_abund_data |>
  dplyr::left_join(abrev, by = "species")
rank_abund_data1 |>
  filter(is.na(abreviacion))
rank_abund_data1
rank_abund_data1 |>
  dplyr::mutate(full_order = seq(1, length(sites))) |>
  tidyr::separate(species, c("genus", "specie"), sep = " ") |>
  dplyr::mutate(
    genus_id = stringr::str_sub(genus, 1, 4),
    spp_id = stringr::str_sub(specie, 1, 4),
    spp_id = dplyr::if_else(stringr::str_detect(specie, "sp.") == TRUE, "sp", spp_id),
    id = paste0(genus_id, " ", spp_id)
  ) |>
  ggplot(aes(full_order, logabun)) +
  geom_point(aes(color = sites, shape = sites),
    size = 2
  ) +
  geom_line(aes(group = sites, color = sites),
    size = 1, linetype = "dotted"
  ) +
  geom_text(aes(label = abreviacion, color = sites),
    hjust = -.2, size = 5,
    angle = 45
  ) +
  # ggrepel::geom_text_repel(aes(label = abreviacion, color = sites),
  #          hjust = -.2, size = 5,
  #          max.overlaps = 2) +
  scale_x_continuous(
    breaks = seq(1, 75, 5),
    limits = c(0, 75)
  ) +
  scale_y_continuous(limits = c(0, 1.8)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    color = "Bosque",
    shape = "Bosque",
    y = "Log10(Abundancia)",
    title = " Ranking de Abundancia "
  )
ggview(
  dpi = 600,
  height = 7,
  width = 11
)
ggsave("img/rank_abund.jpeg",
  dpi = 600,
  height = 7,
  width = 11
)


# Diversidad beta -----------------------------------------------


beta_div <- vegdist(comm_mur, method = "jaccard")

beta_tibble <- beta_div |>
  as.matrix() |>
  corrr::as_cordf() |>
  corrr::stretch(na.rm = TRUE, remove.dups = TRUE) |>
  rename(jaccard = r)
beta_tibble

model <- hclust(beta_div)
ggdendro::ggdendrogram(model, size = 2) +
  theme_bw() +
  coord_flip() +
  labs(
    x = "",
    y = "",
    title = " Dendograma "
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggview(
  dpi = 600,
  units = "cm",
  height = 7,
  width = 12
)
ggsave("img/dendrograma.jpeg",
       dpi = 600,
       units = "cm",
       height = 7,
       width = 12
)


# waffle chart ------------------------------------------------------------
sp_by_site <- data |>
  mutate(tb = case_when(
    tipo_de_bosque == "Bosque de colina baja" ~ "BCB",
    tipo_de_bosque == "Bosque de montaña basimontano" ~ "BMBM",
    tipo_de_bosque == "Bosque de montaña montano" ~ "BMM",
    tipo_de_bosque == "Bosque de terraza baja" ~ "BTB"
  )) |>
  select(-c(fecha, tipo_de_bosque)) |>
  distinct()
sp_by_site
summari_tab <- sp_by_site |>
  rowwise() |>
  mutate(
    group = sum(group_by(sp_by_site, tb) |>
      distinct(especie) |>
      ungroup() |>
      select(especie) |>
      unlist() %in% especie),
    shared = group == length(unique(sp_by_site$tb))
  )

df_dat <- summari_tab |>
  group_by(tb, group) |>
  summarise(
    sp_by_group = n_distinct(especie),
    .groups = "drop"
  ) |>
  mutate(
    tb = factor(tb,
      levels = c("BTB", "BCB", "BMBM", "BMM")
    ),
    group = factor(group, levels = c("1", "2", "3", "4"))
  )

ggplot(df_dat, aes(fill = group, values = sp_by_group)) +
  geom_waffle_p(color = "white", size = .25, n_rows = 10, flip = TRUE) +
  facet_wrap(~tb, nrow = 4, strip.position = "bottom") +
  scale_x_discrete()+
  scale_y_continuous(labels = function(x) x * 10,
                     expand = c(0,0)) +

  theme_minimal() +

  guides(fill = guide_legend(reverse = TRUE)) +
  ggthemes::scale_fill_tableau(name = "Numero de Sitios") +
  coord_equal() +
  labs(
    title = "Especies Presentes",
    x = "Bosque",
    y = "Especies"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(hjust = .5),
    legend.position = "bottom"
  )

ggview(
  dpi = 600,
  units = "cm",
  height = 12,
  width = 8
)
ggsave("img/presence_by_forest.jpeg",
       dpi = 600,
       units = "cm",
       height = 12,
       width = 8
)

