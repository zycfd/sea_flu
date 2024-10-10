## Note

1. **`epi_data/`**: Epidemiological data were downloaded from FluNet (https://www.who.int/tools/flunet) after removing uninformative or non-sensible data entry. We imputed the unsubtyped or lineage-undetermined samples into specific subtypes or lineages based on the available weekly- and country-specific proportion of subtypes or lineages. To do so, we used a Bayesian model with an uninformative Beta (1, 1) prior to calculate the posterior proportions and 95% uncertainty levels (that's why non-integer cases occur in the table). It was described in the methods section of our prior study (https://doi.org/10.1101/2023.12.20.23300299).
2. **`genetic_data/`**: The sequence header IDs in our three sub-sampled datasets are provided.
3. **`geo_data/`**: Latitude and longitude of the centroid of each sub-location are provided.
4. **`others/`**: Metadata of sub-sampled sequences, etc.
