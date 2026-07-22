WITH RECURSIVE gbk_file AS (
  -- Start with each gbk row and its full path as "rest"
  SELECT
    id,
    path,
    path AS rest
  FROM gbk

  UNION ALL

  -- Peel off everything up to the next '/'
  SELECT
    id,
    path,
    SUBSTR(rest, INSTR(rest, '/') + 1) AS rest
  FROM gbk_file
  WHERE INSTR(rest, '/') > 0
),
gbk_basename AS (
  -- When there is no '/' left, rest is the filename
  SELECT
    id,
    rest AS filename
  FROM gbk_file
  WHERE INSTR(rest, '/') = 0
)
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
  gcf_id,
  product,
  sample
ORDER BY
  product,
  sample,
  gcf_id;
