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
' from res
where visited_ids_count = (select max(visited_ids_count) from res)