select group_concat(mol_block || 
'
>  <matched_ids>
' ||
matched_ids ||
'

>  <matched_ids_count>
' ||
matched_ids_count,
'

$$$$
') || '

$$$$
'from mols 
where rowid in (select min(rowid) from mols where visited_ids_count = (select max(visited_ids_count) from mols) group by id)
