use crate::cli::Cli;

#[derive(Clone, Copy)]
pub struct Filter {
    pub min_mapq: u8,
    pub exclude_flags: u16,
    pub include_flags: u16,
    pub count_splice: bool,
}

impl Filter {
    pub fn from_cli(cli: &Cli) -> Self {
        Self {
            min_mapq: cli.min_mapq,
            exclude_flags: cli.exclude_flags,
            include_flags: cli.include_flags,
            count_splice: cli.count_splice,
        }
    }

    #[inline]
    pub fn keep(&self, flags: u16, mapq: u8) -> bool {
        if flags & self.exclude_flags != 0 {
            return false;
        }
        if flags & self.include_flags != self.include_flags {
            return false;
        }
        mapq >= self.min_mapq
    }
}
