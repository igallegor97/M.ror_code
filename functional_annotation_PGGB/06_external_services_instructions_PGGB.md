# Servicios externos: instrucciones y trazabilidad

Use siempre los portales oficiales o Galaxy Europe, anote fecha, versión visible,
parámetros y conserve el archivo original descargado. Los límites de carga pueden
cambiar; divida el FASTA por muestra o en lotes sin renombrar los IDs.

## DeepTMHMM

Para cada muestra, cargue `<muestra>.all_proteins.faa` desde `results/06_external`.
Seleccione predicción de topología para proteínas; descargue el TSV/GFF (y el
archivo de topologías si existe). Guárdelo en:

`results/06_external/deeptmhmm/<muestra>/`

Interpretación para secretoma clásico: conservar SignalP positivo y excluir
proteínas con hélices transmembrana posteriores al péptido señal. No convierta
automáticamente todas las proteínas sin TM en secretadas.

## EffectorP 3

Cargue `<muestra>.signalp_secreted.faa`. Si el FASTA aparece con sufijo `.MISSING`,
revise el nombre de la salida madura de SignalP y cópiela con ese nombre. Elija
modo de proteínas fúngicas cuando el portal lo solicite. Descargue la tabla completa
en `results/06_external/effectorp3/<muestra>/`. Conserve probabilidades/clases, no
solo una lista binaria.

## InterProScan selectivo / Galaxy Europe

1. Añada a `priority_gene_ids_PGGB.txt` los IDs de genes candidatos (p. ej. top 100,
   hotspots, secretados o efectores), uno por línea.
2. Ejecute `bash 05_prepare_external_services_PGGB.sh`.
3. Cargue `priority_candidates.faa` en InterPro o Galaxy Europe.
4. En Galaxy, busque InterProScan, seleccione proteína, TSV + GFF3, y las bases
   disponibles que aporten información adicional (por ejemplo PANTHER, SUPERFAMILY,
   PRINTS, PROSITE). Pfam ya se analizó localmente, por lo que repetirla es opcional.
5. Active GO/pathways si la versión lo permite. Descargue el historial o workflow,
   el TSV y el GFF3 en `results/06_external/interproscan/`.

InterProScan es validación selectiva, no requisito para que finalice la tabla maestra.
