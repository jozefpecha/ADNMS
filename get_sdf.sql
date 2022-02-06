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
'from res 
where rowid in (select min(rowid) from res where visited_ids_count = (select max(visited_ids_count) from res) group by mol_id)
