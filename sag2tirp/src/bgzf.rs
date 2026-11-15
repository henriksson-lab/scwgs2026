//! Standard BGZF, with complete TIRP records kept inside individual blocks.
use std::io::{self, Write};

use noodles::bgzf::{self, VirtualPosition};

const RAW_LIMIT: usize = 60 * 1024;

pub(crate) struct Writer<W: Write> {
    inner: bgzf::io::Writer<W>,
    buffered: usize,
}

impl<W: Write> Writer<W> {
    pub fn new(inner: W) -> Self {
        Self {
            inner: bgzf::io::Writer::new(inner),
            buffered: 0,
        }
    }

    pub fn set_cell(&mut self) -> io::Result<()> {
        self.flush()
    }

    pub fn position(&self) -> VirtualPosition {
        self.inner.virtual_position()
    }

    pub fn reserve_record(&mut self, size: usize) -> io::Result<()> {
        if size > RAW_LIMIT {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "record exceeds BGZF payload limit",
            ));
        }
        if self.buffered + size > RAW_LIMIT {
            self.flush()?;
        }
        Ok(())
    }

    pub fn finish(self) -> io::Result<W> {
        self.inner.finish()
    }
}

impl<W: Write> Write for Writer<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let len = self.inner.write(bytes)?;
        self.buffered += len;
        Ok(len)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()?;
        self.buffered = 0;
        Ok(())
    }
}
