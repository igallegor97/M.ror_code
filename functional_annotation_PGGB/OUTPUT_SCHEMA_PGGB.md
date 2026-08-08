# Esquema normalizado

La clave estable es `(sample_id, protein_id)`.

| Campo | Contenido |
|---|---|
| sample_id | muestra del manifiesto |
| protein_id | primer campo del encabezado FASTA |
| eggnog_description | descripción transferida por eggNOG |
| GO_terms | términos GO de eggNOG |
| KEGG_ko | ortologías KEGG de eggNOG |
| pfam_domains | modelos Pfam-A significativos por `--cut_ga` |
| cazy_families | familias CAZy detectadas en la salida dbCAN |
| signalp_prediction | línea de predicción SignalP 5 |

Las columnas externas se incorporarán tras revisar el formato exacto descargado;
no se mezclan silenciosamente versiones/formats de servidores web.

