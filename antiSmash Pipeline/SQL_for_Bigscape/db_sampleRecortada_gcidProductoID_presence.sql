WITH RECURSIVE gbk_file AS (
  SELECT
    id,
    path AS rest
  FROM gbk

  UNION ALL

  SELECT
    id,
    SUBSTR(rest, INSTR(rest, '/') + 1) AS rest
  FROM gbk_file
  WHERE INSTR(rest, '/') > 0
),
gbk_basename AS (
  SELECT
    id,
    rest AS filename
  FROM gbk_file
  WHERE INSTR(rest, '/') = 0
),
base AS (
  SELECT
    bgc_record_family.family_id AS gcf_id,
    bgc_record.product          AS product,
    gbk_basename.filename       AS sample,
    1                           AS presence
  FROM bgc_record
  JOIN gbk_basename
    ON gbk_basename.id = bgc_record.gbk_id
  JOIN bgc_record_family
    ON bgc_record_family.record_id = bgc_record.id
  WHERE bgc_record.product IS NOT NULL
    AND bgc_record.product <> ''
  GROUP BY
    gcf_id, product, sample
),
with_category AS (
  SELECT
    base.sample,
    base.gcf_id,
    COALESCE(
      "BGCs_tax_table"."Category",
      CASE
        WHEN INSTR(base.product, '.') > 0 THEN
          CASE
            WHEN LOWER(base.product) LIKE '%terpene%' THEN 'terpene'
            WHEN LOWER(base.product) LIKE '%nrps%'    THEN 'NRPS'
            WHEN LOWER(base.product) LIKE '%ripp%'    THEN 'RiPP'
            WHEN LOWER(base.product) LIKE '%pks%'     THEN 'PKS'
            ELSE 'other'
          END
        ELSE 'other'
      END
    ) AS category,
    base.presence
  FROM base
  LEFT JOIN "BGCs_tax_table"
    ON "BGCs_tax_table"."Product" = base.product
)
SELECT
  sample,
  category || '_gcf_' || gcf_id AS gcf_id,
  1 AS presence
FROM with_category
GROUP BY
  sample,
  gcf_id
ORDER BY
  sample,
  gcf_id;
