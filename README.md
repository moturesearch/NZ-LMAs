# NZ-LMAs

This repository contains concordance tables ("crosswalks") between New Zealand area unit and labour market area codes.
We derive these tables from commuting data collected as part of the 2001--2013 New Zealand Censuses.
`paper/` contains the source files for the paper describing our derivation method.

## Data dictionary

Variable | Type | Description
--- | --- | ---
`au13` | int | 2013 area unit code
`lma<yy>` | int | LMA in census year `<yy>`
`slma<yy>` | int | Sub-LMA in census year `<yy>`

## License

[MIT](LICENSE)
